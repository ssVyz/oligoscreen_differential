//! Core analysis algorithms for oligo variant detection

use std::collections::{HashMap, HashSet};
use super::iupac::{base_to_bit, sequence_matches_consensus_bytes, IUPAC_FROM_MASK};
use super::types::{AnalysisMethod, Variant, WindowAnalysisResult};

/// Analyze sequences using the specified method
pub fn analyze_sequences(
    sequences: &[&str],
    method: &AnalysisMethod,
    exclude_n: bool,
    coverage_threshold: f64,
) -> WindowAnalysisResult {
    if sequences.is_empty() {
        return WindowAnalysisResult {
            skipped: true,
            skip_reason: Some("No sequences to analyze".to_string()),
            ..Default::default()
        };
    }

    let total = sequences.len();

    let variants = match method {
        AnalysisMethod::NoAmbiguities => find_variants_no_ambiguities(sequences),
        AnalysisMethod::FixedAmbiguities(max_amb) => {
            find_minimum_variants_greedy(sequences, *max_amb as usize, exclude_n)
        }
        AnalysisMethod::Incremental(target_pct, max_amb) => {
            find_incremental_variants(
                sequences,
                *target_pct as f64,
                exclude_n,
                max_amb.map(|n| n as usize),
            )
        }
    };

    // Calculate variants needed for coverage threshold
    let (variants_needed, coverage_at_threshold) =
        calculate_variants_for_threshold(&variants, total, coverage_threshold);

    WindowAnalysisResult {
        variants,
        total_sequences: total,
        sequences_analyzed: total,
        no_match_count: 0,
        variants_for_threshold: variants_needed,
        coverage_at_threshold,
        skipped: false,
        skip_reason: None,
    }
}

/// Find all unique variants without ambiguity codes
fn find_variants_no_ambiguities(sequences: &[&str]) -> Vec<Variant> {
    let mut counts: HashMap<&str, usize> = HashMap::new();

    for &seq in sequences {
        *counts.entry(seq).or_insert(0) += 1;
    }

    let total = sequences.len() as f64;
    let mut variants: Vec<Variant> = counts
        .into_iter()
        .map(|(seq, count)| Variant {
            sequence: seq.to_string(),
            count,
            percentage: (count as f64 / total) * 100.0,
        })
        .collect();

    // Sort by count descending
    variants.sort_by(|a, b| b.count.cmp(&a.count));
    variants
}

/// Find minimum variants using greedy set cover with ambiguity codes
fn find_minimum_variants_greedy(
    sequences: &[&str],
    max_ambiguities: usize,
    exclude_n: bool,
) -> Vec<Variant> {
    if sequences.is_empty() {
        return Vec::new();
    }

    // Count unique sequences
    let mut seq_counts: HashMap<&str, usize> = HashMap::new();
    for &seq in sequences {
        *seq_counts.entry(seq).or_insert(0) += 1;
    }

    let total = sequences.len() as f64;
    let mut uncovered: HashSet<&str> = seq_counts.keys().copied().collect();
    let mut variants = Vec::new();

    while !uncovered.is_empty() {
        let (best_consensus, best_coverage) = find_best_consensus(
            &uncovered,
            &seq_counts,
            max_ambiguities,
            exclude_n,
        );

        if best_coverage.is_empty() {
            // Fallback: use the most frequent uncovered sequence as-is
            let most_freq = uncovered
                .iter()
                .max_by_key(|&&s| seq_counts.get(s).unwrap_or(&0))
                .copied()
                .unwrap();

            let count = *seq_counts.get(most_freq).unwrap_or(&1);
            variants.push(Variant {
                sequence: most_freq.to_string(),
                count,
                percentage: (count as f64 / total) * 100.0,
            });
            uncovered.remove(most_freq);
        } else {
            let count: usize = best_coverage.iter()
                .map(|&s| seq_counts.get(s).unwrap_or(&0))
                .sum();

            variants.push(Variant {
                sequence: best_consensus,
                count,
                percentage: (count as f64 / total) * 100.0,
            });

            for s in best_coverage {
                uncovered.remove(s);
            }
        }
    }

    variants
}

/// Find the best consensus that covers the most sequences within ambiguity limit.
/// Uses bitmask tracking for zero-allocation inner loop.
fn find_best_consensus<'a>(
    uncovered: &HashSet<&'a str>,
    seq_counts: &HashMap<&'a str, usize>,
    max_ambiguities: usize,
    exclude_n: bool,
) -> (String, HashSet<&'a str>) {
    let mut best_consensus = String::new();
    let mut best_coverage: HashSet<&str> = HashSet::new();
    let mut best_score = 0usize;

    let mut uncovered_sorted: Vec<_> = uncovered.iter().copied().collect();
    uncovered_sorted.sort_by_key(|&s| std::cmp::Reverse(seq_counts.get(s).unwrap_or(&0)));

    let seq_len = uncovered_sorted.first().map(|s| s.len()).unwrap_or(0);
    if seq_len == 0 {
        return (best_consensus, best_coverage);
    }

    let mut group_mask: Vec<u8> = vec![0u8; seq_len];

    for &seed_seq in uncovered_sorted.iter().take(50) {
        // Initialize group_mask from seed
        let seed_bytes = seed_seq.as_bytes();
        for pos in 0..seq_len {
            group_mask[pos] = base_to_bit(seed_bytes[pos]);
        }

        // Try adding other sequences incrementally
        for &other_seq in uncovered {
            if other_seq == seed_seq {
                continue;
            }

            let other_bytes = other_seq.as_bytes();
            let mut trial_amb_count = 0usize;
            let mut trial_valid = true;

            for pos in 0..seq_len {
                let m = group_mask[pos] | base_to_bit(other_bytes[pos]);
                if m.count_ones() > 1 {
                    trial_amb_count += 1;
                    if (exclude_n && m == 0b1111) || trial_amb_count > max_ambiguities {
                        trial_valid = false;
                        break;
                    }
                }
            }

            if trial_valid {
                // Accept: update group_mask in-place
                for pos in 0..seq_len {
                    group_mask[pos] |= base_to_bit(other_bytes[pos]);
                }
            }
        }

        let (consensus, _, is_valid) = consensus_from_mask(&group_mask, exclude_n);
        if !is_valid {
            continue;
        }

        // Check actual coverage
        let consensus_bytes = consensus.as_bytes();
        let mut coverage: HashSet<&str> = HashSet::new();
        for &seq in uncovered {
            if sequence_matches_consensus_bytes(seq.as_bytes(), consensus_bytes) {
                coverage.insert(seq);
            }
        }

        let score: usize = coverage.iter()
            .map(|&s| seq_counts.get(s).unwrap_or(&0))
            .sum();

        if score > best_score {
            best_score = score;
            best_consensus = consensus;
            best_coverage = coverage;
        }
    }

    (best_consensus, best_coverage)
}

/// Find variants incrementally, each covering target percentage of remaining
fn find_incremental_variants(
    sequences: &[&str],
    target_percentage: f64,
    exclude_n: bool,
    max_ambiguities: Option<usize>,
) -> Vec<Variant> {
    if sequences.is_empty() {
        return Vec::new();
    }

    let total_original = sequences.len() as f64;
    let mut remaining: Vec<&str> = sequences.to_vec();
    let mut variants = Vec::new();

    while !remaining.is_empty() {
        let remaining_total = remaining.len();
        let target_count = ((target_percentage / 100.0) * remaining_total as f64).ceil() as usize;

        let mut remaining_counts: HashMap<&str, usize> = HashMap::new();
        for &seq in &remaining {
            *remaining_counts.entry(seq).or_insert(0) += 1;
        }

        let unique_remaining: Vec<&str> = remaining_counts.keys().copied().collect();

        let (best_consensus, best_coverage_count) = find_incremental_consensus(
            &unique_remaining,
            &remaining_counts,
            target_count,
            exclude_n,
            max_ambiguities,
        );

        let percentage = (best_coverage_count as f64 / total_original) * 100.0;
        variants.push(Variant {
            sequence: best_consensus.clone(),
            count: best_coverage_count,
            percentage,
        });

        // Remove covered sequences using byte-level matching
        let best_bytes = best_consensus.as_bytes();
        remaining.retain(|&seq| !sequence_matches_consensus_bytes(seq.as_bytes(), best_bytes));
    }

    variants
}

/// Find consensus for incremental method using bitmask tracking.
fn find_incremental_consensus(
    unique_remaining: &[&str],
    remaining_counts: &HashMap<&str, usize>,
    target_count: usize,
    exclude_n: bool,
    max_ambiguities: Option<usize>,
) -> (String, usize) {
    if unique_remaining.is_empty() {
        return (String::new(), 0);
    }

    let seq_len = unique_remaining[0].len();
    let mut best_consensus = String::new();
    let mut best_coverage_count = 0usize;
    let mut found_target = false;

    let max_amb_level = max_ambiguities.unwrap_or(seq_len);

    // Allocated once per call, reused across all seeds and amb_levels
    let mut group_mask: Vec<u8> = vec![0u8; seq_len];

    for amb_level in 0..=max_amb_level {
        if found_target {
            break;
        }

        let mut sorted_remaining: Vec<_> = unique_remaining.to_vec();
        sorted_remaining.sort_by_key(|&s| std::cmp::Reverse(remaining_counts.get(s).unwrap_or(&0)));

        for &seed_seq in sorted_remaining.iter().take(50) {
            // Initialize group_mask from seed
            let seed_bytes = seed_seq.as_bytes();
            for pos in 0..seq_len {
                group_mask[pos] = base_to_bit(seed_bytes[pos]);
            }

            // Try adding other sequences incrementally via bitmask
            for &other_seq in unique_remaining {
                if other_seq == seed_seq {
                    continue;
                }

                let other_bytes = other_seq.as_bytes();
                let mut trial_amb_count = 0usize;
                let mut trial_valid = true;

                for pos in 0..seq_len {
                    let m = group_mask[pos] | base_to_bit(other_bytes[pos]);
                    if m.count_ones() > 1 {
                        trial_amb_count += 1;
                        if (exclude_n && m == 0b1111) || trial_amb_count > amb_level {
                            trial_valid = false;
                            break;
                        }
                    }
                }

                if trial_valid {
                    for pos in 0..seq_len {
                        group_mask[pos] |= base_to_bit(other_bytes[pos]);
                    }
                }
            }

            let (consensus, amb_count, is_valid) = consensus_from_mask(&group_mask, exclude_n);

            if !is_valid || amb_count > amb_level {
                continue;
            }

            // Check actual coverage using byte-level matching
            let consensus_bytes = consensus.as_bytes();
            let mut coverage_count = 0usize;
            for &seq in unique_remaining {
                if sequence_matches_consensus_bytes(seq.as_bytes(), consensus_bytes) {
                    coverage_count += remaining_counts.get(seq).unwrap_or(&0);
                }
            }

            if coverage_count >= target_count {
                if !found_target || coverage_count > best_coverage_count {
                    best_consensus = consensus;
                    best_coverage_count = coverage_count;
                    found_target = true;
                }
            } else if coverage_count > best_coverage_count {
                best_consensus = consensus;
                best_coverage_count = coverage_count;
            }
        }
    }

    // Fallback
    if best_consensus.is_empty() && !unique_remaining.is_empty() {
        let most_freq = unique_remaining
            .iter()
            .max_by_key(|&&s| remaining_counts.get(s).unwrap_or(&0))
            .copied()
            .unwrap();
        best_consensus = most_freq.to_string();
        best_coverage_count = *remaining_counts.get(most_freq).unwrap_or(&1);
    }

    (best_consensus, best_coverage_count)
}

/// Build a consensus String from a bitmask array.
/// Returns (consensus, ambiguity_count, is_valid).
fn consensus_from_mask(mask: &[u8], exclude_n: bool) -> (String, usize, bool) {
    let mut consensus = String::with_capacity(mask.len());
    let mut amb_count = 0;
    for &m in mask {
        let code = IUPAC_FROM_MASK[m as usize];
        if m.count_ones() > 1 {
            amb_count += 1;
            if exclude_n && code == b'N' {
                return (consensus, amb_count, false);
            }
        }
        consensus.push(code as char);
    }
    (consensus, amb_count, true)
}

/// Create consensus from sequences using bitmask arithmetic.
fn create_consensus_from_seqs(sequences: &[&str], exclude_n: bool) -> (String, usize, bool) {
    if sequences.is_empty() {
        return (String::new(), 0, true);
    }

    let seq_len = sequences[0].len();
    let mut consensus = String::with_capacity(seq_len);
    let mut ambiguity_count = 0;

    for pos in 0..seq_len {
        let mut mask: u8 = 0;
        for &seq in sequences {
            let bytes = seq.as_bytes();
            if pos < bytes.len() {
                mask |= base_to_bit(bytes[pos]);
            }
        }

        let code = IUPAC_FROM_MASK[mask as usize];
        if mask.count_ones() > 1 {
            if exclude_n && code == b'N' {
                return (consensus, ambiguity_count, false);
            }
            ambiguity_count += 1;
        }
        consensus.push(code as char);
    }

    (consensus, ambiguity_count, true)
}

/// Calculate how many variants are needed to reach coverage threshold
fn calculate_variants_for_threshold(
    variants: &[Variant],
    total: usize,
    threshold: f64,
) -> (usize, f64) {
    if variants.is_empty() || total == 0 {
        return (0, 0.0);
    }

    let mut cumulative = 0.0;
    for (i, variant) in variants.iter().enumerate() {
        cumulative += variant.percentage;
        if cumulative >= threshold {
            return (i + 1, cumulative);
        }
    }

    (variants.len(), cumulative)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_no_ambiguities() {
        let seqs = vec!["ACGT", "ACGT", "ACGA", "ACGA", "ACGA"];
        let variants = find_variants_no_ambiguities(&seqs);
        assert_eq!(variants.len(), 2);
        assert_eq!(variants[0].sequence, "ACGA");
        assert_eq!(variants[0].count, 3);
    }

    #[test]
    fn test_calculate_threshold() {
        let variants = vec![
            Variant { sequence: "A".to_string(), count: 50, percentage: 50.0 },
            Variant { sequence: "B".to_string(), count: 30, percentage: 30.0 },
            Variant { sequence: "C".to_string(), count: 20, percentage: 20.0 },
        ];
        let (n, cov) = calculate_variants_for_threshold(&variants, 100, 80.0);
        assert_eq!(n, 2);
        assert_eq!(cov, 80.0);
    }

    #[test]
    fn test_incremental_variants() {
        let seqs = vec!["ACGT", "ACGT", "ACGA", "ACGA", "ACGA", "TCGT", "TCGT"];
        let variants = find_incremental_variants(&seqs, 50.0, false, Some(1));
        assert!(!variants.is_empty());
        let total_count: usize = variants.iter().map(|v| v.count).sum();
        assert_eq!(total_count, 7);
    }

    #[test]
    fn test_fixed_ambiguities() {
        let seqs = vec!["ACGT", "ACGA"];
        let variants = find_minimum_variants_greedy(&seqs, 1, false);
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].count, 2);
    }
}
