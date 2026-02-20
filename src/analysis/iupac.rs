//! IUPAC ambiguity codes and DNA sequence utilities

use std::collections::{HashMap, HashSet};
use once_cell::sync::Lazy;

/// Standard DNA bases
pub const STANDARD_BASES: [char; 4] = ['A', 'C', 'G', 'T'];

/// Ambiguous IUPAC bases (excluding N if needed)
pub const AMBIGUOUS_BASES: [char; 11] = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N'];

/// Gap characters
pub const GAP_CHARS: [char; 2] = ['-', '.'];

/// IUPAC code to bases mapping
pub static IUPAC_TO_BASES: Lazy<HashMap<char, HashSet<char>>> = Lazy::new(|| {
    let mut map = HashMap::new();
    map.insert('A', ['A'].into_iter().collect());
    map.insert('C', ['C'].into_iter().collect());
    map.insert('G', ['G'].into_iter().collect());
    map.insert('T', ['T'].into_iter().collect());
    map.insert('R', ['A', 'G'].into_iter().collect());
    map.insert('Y', ['C', 'T'].into_iter().collect());
    map.insert('S', ['G', 'C'].into_iter().collect());
    map.insert('W', ['A', 'T'].into_iter().collect());
    map.insert('K', ['G', 'T'].into_iter().collect());
    map.insert('M', ['A', 'C'].into_iter().collect());
    map.insert('B', ['C', 'G', 'T'].into_iter().collect());
    map.insert('D', ['A', 'G', 'T'].into_iter().collect());
    map.insert('H', ['A', 'C', 'T'].into_iter().collect());
    map.insert('V', ['A', 'C', 'G'].into_iter().collect());
    map.insert('N', ['A', 'C', 'G', 'T'].into_iter().collect());
    map
});

/// Bases to IUPAC code mapping
pub static BASES_TO_IUPAC: Lazy<HashMap<Vec<char>, char>> = Lazy::new(|| {
    let mut map = HashMap::new();
    map.insert(vec!['A'], 'A');
    map.insert(vec!['C'], 'C');
    map.insert(vec!['G'], 'G');
    map.insert(vec!['T'], 'T');
    map.insert(vec!['A', 'G'], 'R');
    map.insert(vec!['C', 'T'], 'Y');
    map.insert(vec!['C', 'G'], 'S');
    map.insert(vec!['A', 'T'], 'W');
    map.insert(vec!['G', 'T'], 'K');
    map.insert(vec!['A', 'C'], 'M');
    map.insert(vec!['C', 'G', 'T'], 'B');
    map.insert(vec!['A', 'G', 'T'], 'D');
    map.insert(vec!['A', 'C', 'T'], 'H');
    map.insert(vec!['A', 'C', 'G'], 'V');
    map.insert(vec!['A', 'C', 'G', 'T'], 'N');
    map
});

/// Complement mapping for reverse complement
pub static COMPLEMENT: Lazy<HashMap<char, char>> = Lazy::new(|| {
    let mut map = HashMap::new();
    map.insert('A', 'T');
    map.insert('T', 'A');
    map.insert('C', 'G');
    map.insert('G', 'C');
    map.insert('R', 'Y');
    map.insert('Y', 'R');
    map.insert('S', 'S');
    map.insert('W', 'W');
    map.insert('K', 'M');
    map.insert('M', 'K');
    map.insert('B', 'V');
    map.insert('V', 'B');
    map.insert('D', 'H');
    map.insert('H', 'D');
    map.insert('N', 'N');
    map
});

/// Check if a character is a standard DNA base
pub fn is_standard_base(c: char) -> bool {
    matches!(c, 'A' | 'C' | 'G' | 'T')
}

/// Check if a character is an ambiguous base
pub fn is_ambiguous_base(c: char) -> bool {
    matches!(c, 'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V' | 'N')
}

/// Check if a character is a gap
pub fn is_gap(c: char) -> bool {
    matches!(c, '-' | '.')
}

/// Check if a character is a valid DNA character (including ambiguous)
pub fn is_valid_dna(c: char) -> bool {
    is_standard_base(c) || is_ambiguous_base(c)
}

/// Get the IUPAC code for a set of bases
pub fn bases_to_iupac(bases: &HashSet<char>) -> char {
    let mut sorted: Vec<char> = bases.iter().copied().collect();
    sorted.sort();
    *BASES_TO_IUPAC.get(&sorted).unwrap_or(&'N')
}

/// Get the bases represented by an IUPAC code
pub fn iupac_to_bases(code: char) -> Option<&'static HashSet<char>> {
    IUPAC_TO_BASES.get(&code)
}

/// Check if a sequence matches a consensus (with ambiguity codes)
pub fn sequence_matches_consensus(seq: &str, consensus: &str) -> bool {
    if seq.len() != consensus.len() {
        return false;
    }

    for (s, c) in seq.chars().zip(consensus.chars()) {
        if let Some(allowed) = IUPAC_TO_BASES.get(&c) {
            if !allowed.contains(&s) {
                return false;
            }
        } else if s != c {
            return false;
        }
    }
    true
}

/// Create a consensus sequence from multiple sequences
/// Returns (consensus, ambiguity_count, is_valid)
/// is_valid is false if exclude_n=true and N would be required
pub fn create_consensus(sequences: &[&str], exclude_n: bool) -> (String, usize, bool) {
    if sequences.is_empty() {
        return (String::new(), 0, true);
    }

    let seq_len = sequences[0].len();
    let mut consensus = String::with_capacity(seq_len);
    let mut ambiguity_count = 0;

    for pos in 0..seq_len {
        let mut bases_at_pos: HashSet<char> = HashSet::new();
        for seq in sequences {
            if let Some(c) = seq.chars().nth(pos) {
                bases_at_pos.insert(c);
            }
        }

        if bases_at_pos.len() == 1 {
            consensus.push(*bases_at_pos.iter().next().unwrap());
        } else {
            let code = bases_to_iupac(&bases_at_pos);
            if exclude_n && code == 'N' {
                return (consensus, ambiguity_count, false);
            }
            consensus.push(code);
            ambiguity_count += 1;
        }
    }

    (consensus, ambiguity_count, true)
}

/// Compute the reverse complement of a DNA sequence
pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| *COMPLEMENT.get(&c).unwrap_or(&c))
        .collect()
}

/// Count ambiguities in a sequence
pub fn count_ambiguities(seq: &str) -> usize {
    seq.chars().filter(|&c| is_ambiguous_base(c)).count()
}

// ── Bitmask-based IUPAC operations (zero heap allocation) ──────────────────

/// Bitmask representation: bit 0 = A, bit 1 = C, bit 2 = G, bit 3 = T

/// Lookup table: 4-bit bitmask index -> IUPAC code byte.
/// Index 0 (no bases) maps to b'?' and should not occur with valid DNA data.
pub const IUPAC_FROM_MASK: [u8; 16] = [
    b'?', // 0b0000 - no bases (invalid)
    b'A', // 0b0001
    b'C', // 0b0010
    b'M', // 0b0011 - A|C
    b'G', // 0b0100
    b'R', // 0b0101 - A|G
    b'S', // 0b0110 - C|G
    b'V', // 0b0111 - A|C|G
    b'T', // 0b1000
    b'W', // 0b1001 - A|T
    b'Y', // 0b1010 - C|T
    b'H', // 0b1011 - A|C|T
    b'K', // 0b1100 - G|T
    b'D', // 0b1101 - A|G|T
    b'B', // 0b1110 - C|G|T
    b'N', // 0b1111 - A|C|G|T
];

/// Convert a DNA base byte to its bitmask. Also handles IUPAC ambiguity codes.
/// Returns 0 for unrecognized bytes.
#[inline]
pub fn base_to_bit(b: u8) -> u8 {
    match b {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' => 0b1000,
        b'R' => 0b0101,
        b'Y' => 0b1010,
        b'S' => 0b0110,
        b'W' => 0b1001,
        b'K' => 0b1100,
        b'M' => 0b0011,
        b'B' => 0b1110,
        b'D' => 0b1101,
        b'H' => 0b1011,
        b'V' => 0b0111,
        b'N' => 0b1111,
        _ => 0,
    }
}

/// Convert an IUPAC code byte to a bitmask of the bases it represents.
/// Returns 0 for unrecognized bytes.
#[inline]
pub fn iupac_to_mask(b: u8) -> u8 {
    base_to_bit(b)
}

/// Check if a sequence matches a consensus using byte-level bitmask comparison.
/// Zero-allocation equivalent of `sequence_matches_consensus`.
#[inline]
pub fn sequence_matches_consensus_bytes(seq: &[u8], consensus: &[u8]) -> bool {
    if seq.len() != consensus.len() {
        return false;
    }
    for i in 0..seq.len() {
        let base_mask = base_to_bit(seq[i]);
        let cons_mask = iupac_to_mask(consensus[i]);
        if base_mask & cons_mask == 0 {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bitmask_roundtrip() {
        let codes = b"ACGTRYSWKMBDHVN";
        for &code in codes {
            let mask = iupac_to_mask(code);
            assert_eq!(
                IUPAC_FROM_MASK[mask as usize], code,
                "Roundtrip failed for '{}'", code as char
            );
        }
    }

    #[test]
    fn test_base_to_bit() {
        assert_eq!(base_to_bit(b'A'), 0b0001);
        assert_eq!(base_to_bit(b'C'), 0b0010);
        assert_eq!(base_to_bit(b'G'), 0b0100);
        assert_eq!(base_to_bit(b'T'), 0b1000);
        assert_eq!(base_to_bit(b'X'), 0);
    }

    #[test]
    fn test_sequence_matches_consensus_bytes() {
        assert!(sequence_matches_consensus_bytes(b"ACGT", b"ACGT"));
        assert!(sequence_matches_consensus_bytes(b"ACGT", b"NCGT"));
        assert!(sequence_matches_consensus_bytes(b"ACGT", b"RCGT"));
        assert!(!sequence_matches_consensus_bytes(b"ACGT", b"YCGT"));
        assert!(!sequence_matches_consensus_bytes(b"ACG", b"ACGT"));
    }

    #[test]
    fn test_bitmask_matches_hashset_impl() {
        let cases = vec![
            ("ACGT", "ACGT", true),
            ("ACGT", "NCGT", true),
            ("ACGT", "RCGT", true),
            ("TCGT", "RCGT", false),
            ("ACGT", "MCGT", true),
            ("GCGT", "MCGT", false),
        ];
        for (seq, cons, expected) in cases {
            assert_eq!(
                sequence_matches_consensus(seq, cons), expected,
                "HashSet impl: seq={} cons={}", seq, cons
            );
            assert_eq!(
                sequence_matches_consensus_bytes(seq.as_bytes(), cons.as_bytes()), expected,
                "Bitmask impl: seq={} cons={}", seq, cons
            );
        }
    }
}
