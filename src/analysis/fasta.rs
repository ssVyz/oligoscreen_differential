//! FASTA file parsing for template and reference sequences

use super::iupac::{is_ambiguous_base, is_gap, is_standard_base};

/// Parsed template sequence (single sequence)
#[derive(Debug, Clone)]
pub struct TemplateData {
    pub name: String,
    pub sequence: String,
}

/// Parsed reference sequences (multiple, unaligned)
#[derive(Debug, Clone)]
pub struct ReferenceData {
    pub sequences: Vec<String>,
    pub names: Vec<String>,
}

impl ReferenceData {
    pub fn new() -> Self {
        Self {
            sequences: Vec::new(),
            names: Vec::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    pub fn len(&self) -> usize {
        self.sequences.len()
    }
}

impl Default for ReferenceData {
    fn default() -> Self {
        Self::new()
    }
}

/// Parse a single-sequence FASTA as template.
/// Returns error if input contains 0 or more than 1 sequence.
pub fn parse_template_fasta(text: &str) -> Result<TemplateData, String> {
    let (names, sequences) = parse_fasta_sequences(text)?;

    if sequences.is_empty() {
        return Err("No valid sequence found in template input".to_string());
    }
    if sequences.len() > 1 {
        return Err(format!(
            "Template must contain exactly 1 sequence, found {}",
            sequences.len()
        ));
    }

    // Validate template has only standard bases (no gaps or ambiguities)
    let seq = &sequences[0];
    for (i, c) in seq.chars().enumerate() {
        if !is_standard_base(c) {
            return Err(format!(
                "Template contains invalid character '{}' at position {}. Only A, C, G, T are allowed.",
                c, i + 1
            ));
        }
    }

    Ok(TemplateData {
        name: names[0].clone(),
        sequence: sequences[0].clone(),
    })
}

/// Parse multi-sequence FASTA as reference set (unaligned, no length normalization).
pub fn parse_reference_fasta(text: &str) -> Result<ReferenceData, String> {
    let (names, sequences) = parse_fasta_sequences(text)?;

    if sequences.is_empty() {
        return Err("No valid sequences found in reference input".to_string());
    }

    let mut data = ReferenceData::new();
    data.names = names;
    data.sequences = sequences;
    Ok(data)
}

/// Core FASTA parsing: extract names and sequences from FASTA text.
/// Does NOT normalize lengths (suitable for unaligned sequences).
fn parse_fasta_sequences(text: &str) -> Result<(Vec<String>, Vec<String>), String> {
    let mut names = Vec::new();
    let mut sequences = Vec::new();
    let mut current_name = String::new();
    let mut current_seq = String::new();

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        if let Some(name) = line.strip_prefix('>') {
            // Save previous sequence if exists
            if !current_seq.is_empty() {
                names.push(current_name.clone());
                sequences.push(current_seq.clone());
                current_seq.clear();
            }
            current_name = name.to_string();
        } else {
            // Append to current sequence, converting to uppercase
            for c in line.chars() {
                let c = c.to_ascii_uppercase();
                if is_standard_base(c) || is_ambiguous_base(c) || is_gap(c) {
                    if c == '.' {
                        current_seq.push('-');
                    } else {
                        current_seq.push(c);
                    }
                }
                // Ignore other characters (whitespace, numbers, etc.)
            }
        }
    }

    // Don't forget the last sequence
    if !current_seq.is_empty() {
        if current_name.is_empty() {
            current_name = format!("Sequence_{}", sequences.len() + 1);
        }
        names.push(current_name);
        sequences.push(current_seq);
    }

    // If no FASTA headers found, try treating each line as a sequence
    if sequences.is_empty() {
        for (i, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('>') {
                continue;
            }

            let mut seq = String::new();
            for c in line.chars() {
                let c = c.to_ascii_uppercase();
                if is_standard_base(c) || is_ambiguous_base(c) || is_gap(c) {
                    if c == '.' {
                        seq.push('-');
                    } else {
                        seq.push(c);
                    }
                }
            }

            if !seq.is_empty() {
                names.push(format!("Sequence_{}", i + 1));
                sequences.push(seq);
            }
        }
    }

    Ok((names, sequences))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_template() {
        let fasta = ">Template\nACGTACGT";
        let data = parse_template_fasta(fasta).unwrap();
        assert_eq!(data.name, "Template");
        assert_eq!(data.sequence, "ACGTACGT");
    }

    #[test]
    fn test_parse_template_rejects_multiple() {
        let fasta = ">Seq1\nACGT\n>Seq2\nACGT";
        assert!(parse_template_fasta(fasta).is_err());
    }

    #[test]
    fn test_parse_template_rejects_gaps() {
        let fasta = ">Template\nAC-TACGT";
        assert!(parse_template_fasta(fasta).is_err());
    }

    #[test]
    fn test_parse_references() {
        let fasta = ">Ref1\nACGTACGT\n>Ref2\nACGTACGTTT\n>Ref3\nACGT";
        let data = parse_reference_fasta(fasta).unwrap();
        assert_eq!(data.len(), 3);
        // Sequences should NOT be padded to same length
        assert_eq!(data.sequences[0].len(), 8);
        assert_eq!(data.sequences[1].len(), 10);
        assert_eq!(data.sequences[2].len(), 4);
    }
}
