//! Data types for oligo analysis

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Analysis method selection
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnalysisMethod {
    /// Find all unique variants without using ambiguity codes
    NoAmbiguities,
    /// Find minimum variants using up to N ambiguity codes per variant
    FixedAmbiguities(u32),
    /// Incremental: find variants covering X% of remaining sequences each step
    /// Parameters: (target_percentage, optional_max_ambiguities)
    Incremental(u32, Option<u32>),
}

impl Default for AnalysisMethod {
    fn default() -> Self {
        Self::NoAmbiguities
    }
}

impl AnalysisMethod {
    pub fn description(&self) -> String {
        match self {
            Self::NoAmbiguities => "No Ambiguities (exact variants only)".to_string(),
            Self::FixedAmbiguities(n) => format!("Fixed Ambiguities (max {} per variant)", n),
            Self::Incremental(pct, _) => format!("Incremental ({}% coverage per step)", pct),
        }
    }
}

/// Thread count configuration
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ThreadCount {
    /// Use all available CPU cores
    Auto,
    /// Use a specific number of threads
    Fixed(usize),
}

impl Default for ThreadCount {
    fn default() -> Self {
        Self::Auto
    }
}

impl ThreadCount {
    /// Get the actual number of threads to use
    pub fn get_count(&self) -> usize {
        match self {
            Self::Auto => std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1),
            Self::Fixed(n) => *n,
        }
    }
}

/// Pairwise alignment parameters
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PairwiseParams {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open_penalty: i32,
    pub gap_extend_penalty: i32,
    pub max_mismatches: u32,
}

impl Default for PairwiseParams {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -1,
            gap_open_penalty: -2,
            gap_extend_penalty: -1,
            max_mismatches: 8,
        }
    }
}

/// Global analysis parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisParams {
    pub method: AnalysisMethod,
    pub pairwise: PairwiseParams,
    pub exclude_n: bool,
    pub min_oligo_length: u32,
    pub max_oligo_length: u32,
    pub resolution: u32,
    pub coverage_threshold: f64,
    pub thread_count: ThreadCount,
}

impl Default for AnalysisParams {
    fn default() -> Self {
        Self {
            method: AnalysisMethod::NoAmbiguities,
            pairwise: PairwiseParams::default(),
            exclude_n: true,
            min_oligo_length: 18,
            max_oligo_length: 25,
            resolution: 1,
            coverage_threshold: 90.0,
            thread_count: ThreadCount::Auto,
        }
    }
}

/// A single variant with its count and percentage
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Variant {
    pub sequence: String,
    pub count: usize,
    pub percentage: f64,
}

/// Result of analyzing a single window position
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WindowAnalysisResult {
    pub variants: Vec<Variant>,
    pub total_sequences: usize,
    pub sequences_analyzed: usize,
    pub no_match_count: usize,
    pub variants_for_threshold: usize,
    pub coverage_at_threshold: f64,
    pub skipped: bool,
    pub skip_reason: Option<String>,
}

impl Default for WindowAnalysisResult {
    fn default() -> Self {
        Self {
            variants: Vec::new(),
            total_sequences: 0,
            sequences_analyzed: 0,
            no_match_count: 0,
            variants_for_threshold: 0,
            coverage_at_threshold: 0.0,
            skipped: false,
            skip_reason: None,
        }
    }
}

/// Result for a specific oligo length across all positions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LengthResult {
    pub oligo_length: u32,
    pub positions: Vec<PositionResult>,
}

/// Result at a specific template position
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PositionResult {
    pub position: usize,
    pub variants_needed: usize,
    pub analysis: WindowAnalysisResult,
    #[serde(default)]
    pub exclusivity: Option<ExclusivityResult>,
}

/// Exclusivity analysis result for a single position/length
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExclusivityResult {
    pub total_sequences: usize,
    pub no_match_count: usize,
    /// Histogram of mismatch counts, sorted ascending by mismatch count
    pub mismatch_histogram: Vec<MismatchBucket>,
    /// Minimum mismatches across all exclusivity sequences (None = all are no-match)
    pub min_mismatches: Option<u32>,
}

/// A single bucket in the mismatch histogram
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MismatchBucket {
    pub mismatches: u32,
    pub count: usize,
    pub example_name: String,
}

/// Complete screening results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScreeningResults {
    pub params: AnalysisParams,
    pub template_length: usize,
    pub total_sequences: usize,
    pub template_sequence: String,
    pub results_by_length: HashMap<u32, LengthResult>,
    #[serde(default)]
    pub differential_enabled: bool,
    #[serde(default)]
    pub exclusivity_sequence_count: Option<usize>,
}

impl ScreeningResults {
    pub fn new(
        params: AnalysisParams,
        template_length: usize,
        total_sequences: usize,
        template_sequence: String,
        differential_enabled: bool,
        exclusivity_sequence_count: Option<usize>,
    ) -> Self {
        Self {
            params,
            template_length,
            total_sequences,
            template_sequence,
            results_by_length: HashMap::new(),
            differential_enabled,
            exclusivity_sequence_count,
        }
    }
}

/// Progress update during analysis
#[derive(Debug, Clone)]
pub struct ProgressUpdate {
    pub current_length: u32,
    pub current_position: usize,
    pub total_positions: usize,
    pub lengths_completed: u32,
    pub total_lengths: u32,
    pub message: String,
}
