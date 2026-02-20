//! Pairwise alignment logic for matching oligos against reference sequences
//!
//! Uses Smith-Waterman local alignment from the bio crate to find the best
//! match for each template oligo in each reference sequence.

use bio::alignment::pairwise::{Aligner, MatchFunc, MatchParams};
use bio::alignment::AlignmentOperation;

use super::types::PairwiseParams;

/// Concrete Aligner type using MatchParams (nameable, unlike closure-based Aligners).
pub type DnaAligner = Aligner<MatchParams>;

/// Create an Aligner sized for the given dimensions.
pub fn create_aligner(
    oligo_len: usize,
    max_ref_len: usize,
    params: &PairwiseParams,
) -> DnaAligner {
    let match_fn = MatchParams::new(params.match_score, params.mismatch_score);
    Aligner::with_capacity(
        oligo_len,
        max_ref_len,
        params.gap_open_penalty,
        params.gap_extend_penalty,
        match_fn,
    )
}

/// Result of aligning an oligo against a single reference sequence
#[derive(Debug, Clone)]
pub struct PairwiseMatch {
    /// The matched region extracted from the reference (gap-free)
    pub matched_sequence: String,
    /// Alignment score
    pub score: i32,
    /// Number of mismatches in the alignment
    pub mismatches: usize,
    /// Whether the alignment contains gaps (insertions or deletions)
    pub has_gaps: bool,
    /// Whether the alignment covers the full query (oligo)
    pub full_coverage: bool,
}

/// Process an alignment result from a pre-existing aligner.
/// Shared logic used by both standalone alignment and batch collection.
fn process_alignment<F: MatchFunc>(
    aligner: &mut Aligner<F>,
    oligo: &[u8],
    reference: &[u8],
) -> PairwiseMatch {
    let alignment = aligner.local(oligo, reference);

    let mut has_gaps = false;
    let mut mismatches = 0;

    for op in &alignment.operations {
        match op {
            AlignmentOperation::Match => {}
            AlignmentOperation::Subst => {
                mismatches += 1;
            }
            AlignmentOperation::Del | AlignmentOperation::Ins => {
                has_gaps = true;
            }
            AlignmentOperation::Xclip(_) | AlignmentOperation::Yclip(_) => {}
        }
    }

    let aligned_query_len = alignment.xend - alignment.xstart;
    let full_coverage = aligned_query_len == oligo.len();

    let matched_sequence = if !has_gaps && full_coverage {
        String::from_utf8_lossy(&reference[alignment.ystart..alignment.yend]).to_string()
    } else {
        String::new()
    };

    PairwiseMatch {
        matched_sequence,
        score: alignment.score,
        mismatches,
        has_gaps,
        full_coverage,
    }
}

/// Align an oligo against a single reference sequence using local alignment.
/// Creates its own aligner — use `collect_matches` for batch alignment.
pub fn align_oligo_to_reference(
    oligo: &[u8],
    reference: &[u8],
    params: &PairwiseParams,
) -> PairwiseMatch {
    let match_score = params.match_score;
    let mismatch_score = params.mismatch_score;

    let mut aligner = Aligner::with_capacity(
        oligo.len(),
        reference.len(),
        params.gap_open_penalty,
        params.gap_extend_penalty,
        |a: u8, b: u8| -> i32 {
            if a == b { match_score } else { mismatch_score }
        },
    );

    process_alignment(&mut aligner, oligo, reference)
}

/// Align an oligo against all reference sequences and collect valid matches.
///
/// Creates a single aligner sized for the longest reference and reuses it
/// for all alignments, avoiding repeated large allocations.
///
/// Returns (matched_sequences, no_match_count).
/// A match is rejected (counted as "no match") if:
/// - The alignment contains gaps
/// - The alignment doesn't cover the full oligo
/// - The number of mismatches exceeds max_mismatches
pub fn collect_matches(
    oligo: &[u8],
    references: &[Vec<u8>],
    params: &PairwiseParams,
) -> (Vec<String>, usize) {
    let mut matched = Vec::new();
    let mut no_match_count = 0;

    if references.is_empty() {
        return (matched, no_match_count);
    }

    // Create a single aligner sized for the longest reference, reused for all alignments.
    // This avoids re-allocating the O(m*n) DP matrices for every reference.
    let max_ref_len = references.iter().map(|r| r.len()).max().unwrap();
    let match_score = params.match_score;
    let mismatch_score = params.mismatch_score;

    let mut aligner = Aligner::with_capacity(
        oligo.len(),
        max_ref_len,
        params.gap_open_penalty,
        params.gap_extend_penalty,
        |a: u8, b: u8| -> i32 {
            if a == b { match_score } else { mismatch_score }
        },
    );

    for reference in references {
        let result = process_alignment(&mut aligner, oligo, reference);

        if !result.full_coverage || result.has_gaps || result.mismatches > params.max_mismatches as usize
        {
            no_match_count += 1;
        } else {
            matched.push(result.matched_sequence);
        }
    }

    (matched, no_match_count)
}

/// Align an oligo against all references using a pre-existing aligner.
/// The aligner must be sized for at least (oligo.len(), max_ref_len).
pub fn collect_matches_with_aligner(
    aligner: &mut DnaAligner,
    oligo: &[u8],
    references: &[Vec<u8>],
    params: &PairwiseParams,
) -> (Vec<String>, usize) {
    let mut matched = Vec::new();
    let mut no_match_count = 0;

    for reference in references {
        let result = process_alignment(aligner, oligo, reference);

        if !result.full_coverage || result.has_gaps || result.mismatches > params.max_mismatches as usize
        {
            no_match_count += 1;
        } else {
            matched.push(result.matched_sequence);
        }
    }

    (matched, no_match_count)
}

/// Align an oligo against all references using a pre-existing aligner and
/// return per-sequence mismatch counts for exclusivity analysis.
///
/// Returns a Vec with one entry per reference: Some(mismatches) for valid
/// alignments, None for no-match (gaps, partial coverage, or exceeds max_mismatches).
pub fn collect_mismatch_counts_with_aligner(
    aligner: &mut DnaAligner,
    oligo: &[u8],
    references: &[Vec<u8>],
    params: &PairwiseParams,
) -> Vec<Option<u32>> {
    references
        .iter()
        .map(|reference| {
            let result = process_alignment(aligner, oligo, reference);
            if !result.full_coverage
                || result.has_gaps
                || result.mismatches > params.max_mismatches as usize
            {
                None
            } else {
                Some(result.mismatches as u32)
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_params() -> PairwiseParams {
        PairwiseParams::default()
    }

    #[test]
    fn test_exact_match() {
        let oligo = b"TATGGTACGT";
        let reference = b"TATGGTACGTCATGTTCTAGAAATGGGCTGT";
        let result = align_oligo_to_reference(oligo, reference, &default_params());

        assert!(!result.has_gaps);
        assert!(result.full_coverage);
        assert_eq!(result.mismatches, 0);
        assert_eq!(result.matched_sequence, "TATGGTACGT");
    }

    #[test]
    fn test_match_with_offset() {
        let oligo = b"TATGGTACGT";
        let reference = b"AATATGGTACGTCATGTTCTAGAAATGGGCTGT";
        let result = align_oligo_to_reference(oligo, reference, &default_params());

        assert!(!result.has_gaps);
        assert!(result.full_coverage);
        assert_eq!(result.mismatches, 0);
        assert_eq!(result.matched_sequence, "TATGGTACGT");
    }

    #[test]
    fn test_match_with_mismatch() {
        let oligo = b"TATGGTACGT";
        let reference = b"TATGGTTCGTCATGTTCTAGAAATGGGCTGTTTT";
        let result = align_oligo_to_reference(oligo, reference, &default_params());

        assert!(!result.has_gaps);
        assert!(result.full_coverage);
        assert_eq!(result.mismatches, 1);
        assert_eq!(result.matched_sequence, "TATGGTTCGT");
    }

    #[test]
    fn test_collect_matches_from_example() {
        let oligo = b"TATGGTACGT";
        let references: Vec<Vec<u8>> = vec![
            b"TATGGTACGTCATGTTCTAGAAATGGGCTGT".to_vec(),
            b"AATATGGTACGTCATGTTCTAGAAATGGGCTGT".to_vec(),
            b"TATGGTTCGTCATGTTCTAGAAATGGGCTGTTTT".to_vec(),
            b"GTATGGTACGTCATGTTCTAGAAATGGGCTGT".to_vec(),
        ];
        let params = default_params();
        let (matched, no_match) = collect_matches(oligo, &references, &params);

        assert_eq!(no_match, 0); // All should match (1 mismatch is within default max)
        assert_eq!(matched.len(), 4);
        // 3 exact + 1 with mismatch
        assert_eq!(matched.iter().filter(|s| *s == "TATGGTACGT").count(), 3);
        assert_eq!(matched.iter().filter(|s| *s == "TATGGTTCGT").count(), 1);
    }

    #[test]
    fn test_max_mismatches_filter() {
        let oligo = b"TATGGTACGT";
        let references: Vec<Vec<u8>> = vec![
            b"TATGGTACGTCATGTT".to_vec(),
            b"TATGGTTCGTCATGTT".to_vec(), // 1 mismatch
        ];
        let mut params = default_params();
        params.max_mismatches = 0; // No mismatches allowed

        let (matched, no_match) = collect_matches(oligo, &references, &params);
        assert_eq!(matched.len(), 1);
        assert_eq!(no_match, 1);
    }
}
