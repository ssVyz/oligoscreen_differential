//! Main application state and UI

use eframe::egui;
use std::sync::mpsc::{channel, Receiver};
use std::thread;

use crate::analysis::{
    parse_reference_fasta, parse_template_fasta, reverse_complement, run_screening,
    AnalysisMethod, AnalysisParams, ProgressUpdate, ReferenceData, ScreeningResults, TemplateData,
    ThreadCount,
};

/// Info about an imported exclusivity file (UI-only, not serialized)
struct ExclusivityFileEntry {
    file_name: String,
    file_content: String,
    sequence_count: usize,
    min_length: usize,
    max_length: usize,
}

/// Application state
pub struct OligoscreenApp {
    // Input tab state - template
    template_file_name: Option<String>,
    template_data: Option<TemplateData>,
    template_error: Option<String>,

    // Input tab state - references
    reference_file_name: Option<String>,
    reference_data: Option<ReferenceData>,
    reference_error: Option<String>,

    // Differential analysis input
    use_differential: bool,
    exclusivity_files: Vec<ExclusivityFileEntry>,
    exclusivity_data: Option<ReferenceData>,
    exclusivity_error: Option<String>,

    // Analysis parameters
    params: AnalysisParams,
    method_selection: MethodSelection,
    thread_selection: ThreadSelection,
    manual_thread_count: usize,

    // Incremental method options
    incremental_limit_ambiguities: bool,
    incremental_max_ambiguities: u32,

    // Analysis state
    is_analyzing: bool,
    analysis_progress: Option<ProgressUpdate>,
    progress_rx: Option<Receiver<ProgressUpdate>>,
    results_rx: Option<Receiver<ScreeningResults>>,

    // Results state
    results: Option<ScreeningResults>,
    selected_position: Option<usize>,
    selected_length_for_detail: Option<u32>,
    show_detail_window: bool,

    // Detail window display options
    detail_show_reverse_complement: bool,
    detail_show_codon_spacing: bool,

    // View state
    current_tab: Tab,
    zoom_level: f32,

    // Results viewer settings (adjustable without re-running analysis)
    view_coverage_threshold: f64,
    color_green_at: usize,
    color_red_at: usize,
    nomatch_ok_percent: f64,
    nomatch_bad_percent: f64,

    // Differential mode display settings
    differential_mode: bool,
    diff_green_at: u32,
    diff_red_at: u32,
    diff_ignore_count: usize,

    // Save/Load
    save_error: Option<String>,
    load_error: Option<String>,

    // Deferred actions
    pending_save: bool,
    pending_remove_excl: Option<usize>,

    // Output folder for auto-save
    output_folder: Option<String>,

    // Worklist
    next_job_id: u64,
    worklist: Vec<WorklistJob>,
    completed_jobs: Vec<CompletedJob>,
    worklist_state: WorklistState,
    current_job_index: usize,
    selected_completed_job_index: Option<usize>,
    auto_save_error: Option<String>,
    /// Total jobs at the start of a processing batch (for overall progress bar)
    worklist_total_at_start: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Tab {
    Input,
    Analysis,
    Worklist,
    Results,
}

/// A single job in the worklist queue.
/// Captures all inputs and analysis parameters at the time of "Add to Worklist".
struct WorklistJob {
    id: u64,
    // Captured inputs
    template_file_name: String,
    template_data: TemplateData,
    reference_file_name: String,
    reference_data: ReferenceData,
    use_differential: bool,
    exclusivity_file_names: Vec<String>,
    exclusivity_data: Option<ReferenceData>,
    // Captured params (fully resolved method, thread count applied at run time)
    params: AnalysisParams,
    // Output folder (optional, for auto-save)
    output_folder: Option<String>,
    // Summary info for display
    template_length: usize,
    reference_count: usize,
    exclusivity_count: usize,
}

/// A completed job with its results.
struct CompletedJob {
    job: WorklistJob,
    results: ScreeningResults,
}

/// Worklist processing state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WorklistState {
    Idle,
    Processing,
    StopRequested,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MethodSelection {
    NoAmbiguities,
    FixedAmbiguities,
    Incremental,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ThreadSelection {
    Auto,
    Manual,
}

impl Default for OligoscreenApp {
    fn default() -> Self {
        let available_threads = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1);
        Self {
            template_file_name: None,
            template_data: None,
            template_error: None,
            reference_file_name: None,
            reference_data: None,
            reference_error: None,
            use_differential: false,
            exclusivity_files: Vec::new(),
            exclusivity_data: None,
            exclusivity_error: None,
            params: AnalysisParams::default(),
            method_selection: MethodSelection::NoAmbiguities,
            thread_selection: ThreadSelection::Auto,
            manual_thread_count: available_threads,
            incremental_limit_ambiguities: false,
            incremental_max_ambiguities: 3,
            is_analyzing: false,
            analysis_progress: None,
            progress_rx: None,
            results_rx: None,
            results: None,
            selected_position: None,
            selected_length_for_detail: None,
            show_detail_window: false,
            detail_show_reverse_complement: false,
            detail_show_codon_spacing: true,
            current_tab: Tab::Input,
            zoom_level: 1.0,
            view_coverage_threshold: 95.0,
            color_green_at: 1,
            color_red_at: 10,
            nomatch_ok_percent: 5.0,
            nomatch_bad_percent: 50.0,
            differential_mode: false,
            diff_green_at: 5,
            diff_red_at: 0,
            diff_ignore_count: 0,
            save_error: None,
            load_error: None,
            pending_save: false,
            pending_remove_excl: None,
            output_folder: None,
            next_job_id: 1,
            worklist: Vec::new(),
            completed_jobs: Vec::new(),
            worklist_state: WorklistState::Idle,
            current_job_index: 0,
            selected_completed_job_index: None,
            auto_save_error: None,
            worklist_total_at_start: 0,
        }
    }
}

impl OligoscreenApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self::default()
    }

    /// Recalculate variants_for_threshold and coverage_at_threshold for all
    /// positions using the current view_coverage_threshold, without re-running
    /// the full analysis.
    fn recalculate_coverage_threshold(&mut self) {
        let threshold = self.view_coverage_threshold;
        let Some(results) = &mut self.results else {
            return;
        };

        for length_result in results.results_by_length.values_mut() {
            for pos_result in &mut length_result.positions {
                if pos_result.analysis.skipped {
                    continue;
                }
                let mut cumulative = 0.0;
                let mut new_needed = pos_result.analysis.variants.len();
                let mut new_coverage = 0.0;
                for (i, variant) in pos_result.analysis.variants.iter().enumerate() {
                    cumulative += variant.percentage;
                    if cumulative >= threshold {
                        new_needed = i + 1;
                        new_coverage = cumulative;
                        break;
                    }
                }
                if cumulative < threshold {
                    new_coverage = cumulative;
                }
                pos_result.analysis.variants_for_threshold = new_needed;
                pos_result.analysis.coverage_at_threshold = new_coverage;
                pos_result.variants_needed = new_needed;
            }
        }
    }

    /// Resolve the current UI method selection into a concrete AnalysisMethod.
    fn resolve_method(&self) -> AnalysisMethod {
        match self.method_selection {
            MethodSelection::NoAmbiguities => AnalysisMethod::NoAmbiguities,
            MethodSelection::FixedAmbiguities => {
                AnalysisMethod::FixedAmbiguities(self.params.method.get_fixed_ambiguities())
            }
            MethodSelection::Incremental => {
                let max_amb = if self.incremental_limit_ambiguities {
                    Some(self.incremental_max_ambiguities)
                } else {
                    None
                };
                AnalysisMethod::Incremental(self.params.method.get_incremental_pct(), max_amb)
            }
        }
    }

    /// Capture current inputs + params into a WorklistJob and clear the inputs.
    fn add_to_worklist(&mut self) {
        let Some(template_data) = self.template_data.clone() else {
            return;
        };
        let Some(reference_data) = self.reference_data.clone() else {
            return;
        };

        let template_file_name = self.template_file_name.clone().unwrap_or_default();
        let reference_file_name = self.reference_file_name.clone().unwrap_or_default();

        let mut params = self.params.clone();
        params.method = self.resolve_method();

        let exclusivity_file_names: Vec<String> = self
            .exclusivity_files
            .iter()
            .map(|e| e.file_name.clone())
            .collect();
        let exclusivity_data = if self.use_differential {
            self.exclusivity_data.clone()
        } else {
            None
        };

        let template_length = template_data.sequence.len();
        let reference_count = reference_data.len();
        let exclusivity_count = exclusivity_data.as_ref().map(|d| d.len()).unwrap_or(0);

        let job = WorklistJob {
            id: self.next_job_id,
            template_file_name,
            template_data,
            reference_file_name,
            reference_data,
            use_differential: self.use_differential,
            exclusivity_file_names,
            exclusivity_data,
            params,
            output_folder: self.output_folder.clone(),
            template_length,
            reference_count,
            exclusivity_count,
        };

        self.next_job_id += 1;
        self.worklist.push(job);

        // Clear input fields for next job
        self.template_file_name = None;
        self.template_data = None;
        self.template_error = None;
        self.reference_file_name = None;
        self.reference_data = None;
        self.reference_error = None;
        self.exclusivity_files.clear();
        self.exclusivity_data = None;
        self.exclusivity_error = None;
        self.use_differential = false;
    }

    fn select_output_folder(&mut self) {
        if let Some(path) = rfd::FileDialog::new().pick_folder() {
            self.output_folder = Some(path.to_string_lossy().to_string());
        }
    }

    fn remove_worklist_job(&mut self, index: usize) {
        if index < self.worklist.len() {
            // Don't allow removing the currently-processing job
            if self.worklist_state == WorklistState::Processing && index == self.current_job_index {
                return;
            }
            self.worklist.remove(index);
            if self.worklist_state == WorklistState::Processing && index < self.current_job_index {
                self.current_job_index -= 1;
            }
        }
    }

    fn start_worklist_processing(&mut self) {
        if self.worklist.is_empty() || self.worklist_state == WorklistState::Processing {
            return;
        }
        self.worklist_state = WorklistState::Processing;
        self.current_job_index = 0;
        self.worklist_total_at_start = self.worklist.len();
        self.start_next_job();
    }

    fn start_next_job(&mut self) {
        if self.current_job_index >= self.worklist.len() {
            self.worklist_state = WorklistState::Idle;
            self.analysis_progress = None;
            return;
        }

        if self.worklist_state == WorklistState::StopRequested {
            self.worklist_state = WorklistState::Idle;
            self.analysis_progress = None;
            return;
        }

        let job = &self.worklist[self.current_job_index];

        // Apply thread count from Worklist tab controls (not from job snapshot)
        let mut params = job.params.clone();
        params.thread_count = match self.thread_selection {
            ThreadSelection::Auto => ThreadCount::Auto,
            ThreadSelection::Manual => ThreadCount::Fixed(self.manual_thread_count),
        };

        let template_clone = job.template_data.clone();
        let references_clone = job.reference_data.clone();
        let exclusivity_clone = job.exclusivity_data.clone();

        let (progress_tx, progress_rx) = channel();
        let (results_tx, results_rx) = channel();

        self.progress_rx = Some(progress_rx);
        self.results_rx = Some(results_rx);
        self.is_analyzing = true;
        self.analysis_progress = None;

        thread::spawn(move || {
            let results = run_screening(
                &template_clone,
                &references_clone,
                &params,
                exclusivity_clone.as_ref(),
                Some(progress_tx),
            );
            let _ = results_tx.send(results);
        });
    }

    fn check_analysis_progress(&mut self) {
        if let Some(rx) = &self.progress_rx {
            while let Ok(progress) = rx.try_recv() {
                self.analysis_progress = Some(progress);
            }
        }

        if let Some(rx) = &self.results_rx {
            if let Ok(results) = rx.try_recv() {
                self.is_analyzing = false;
                self.progress_rx = None;
                self.results_rx = None;

                // Remove the completed job from the worklist
                let job = self.worklist.remove(self.current_job_index);

                // Auto-save if output folder is set
                if let Some(ref folder) = job.output_folder {
                    let folder = folder.clone();
                    self.auto_save_results(&results, &folder, &job);
                }

                self.completed_jobs.push(CompletedJob { job, results });

                // Select the newly completed job for viewing
                let idx = self.completed_jobs.len() - 1;
                self.selected_completed_job_index = Some(idx);
                self.results = Some(self.completed_jobs[idx].results.clone());
                self.view_coverage_threshold =
                    self.completed_jobs[idx].results.params.coverage_threshold;
                self.differential_mode = self.completed_jobs[idx].results.differential_enabled;

                // current_job_index stays the same because we removed the element at it
                self.start_next_job();
            }
        }
    }

    fn auto_save_results(
        &mut self,
        results: &ScreeningResults,
        folder: &str,
        job: &WorklistJob,
    ) {
        let sanitized_name: String = job
            .template_file_name
            .chars()
            .map(|c| {
                if c.is_alphanumeric() || c == '-' || c == '_' || c == '.' {
                    c
                } else {
                    '_'
                }
            })
            .collect();
        let file_name = format!("{}_{}.json", sanitized_name, job.id);
        let path = std::path::Path::new(folder).join(file_name);

        match serde_json::to_string_pretty(results) {
            Ok(json) => {
                if let Err(e) = std::fs::write(&path, json) {
                    self.auto_save_error = Some(format!("Auto-save failed: {}", e));
                } else {
                    self.auto_save_error = None;
                }
            }
            Err(e) => {
                self.auto_save_error = Some(format!("Auto-save serialize failed: {}", e));
            }
        }
    }

    fn save_results(&mut self) {
        let Some(results) = &self.results else {
            self.save_error = Some("No results to save".to_string());
            return;
        };

        if let Some(path) = rfd::FileDialog::new()
            .add_filter("JSON", &["json"])
            .set_file_name("screening_results.json")
            .save_file()
        {
            match serde_json::to_string_pretty(results) {
                Ok(json) => {
                    if let Err(e) = std::fs::write(&path, json) {
                        self.save_error = Some(format!("Failed to write file: {}", e));
                    } else {
                        self.save_error = None;
                    }
                }
                Err(e) => {
                    self.save_error = Some(format!("Failed to serialize: {}", e));
                }
            }
        }
    }

    fn load_results_into_completed(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("JSON", &["json"])
            .pick_file()
        {
            match std::fs::read_to_string(&path) {
                Ok(json) => match serde_json::from_str::<ScreeningResults>(&json) {
                    Ok(results) => {
                        let file_name = path
                            .file_name()
                            .map(|n| n.to_string_lossy().to_string())
                            .unwrap_or_else(|| "loaded".to_string());

                        let job = WorklistJob {
                            id: self.next_job_id,
                            template_file_name: format!("(loaded) {}", file_name),
                            template_data: TemplateData {
                                name: "Loaded".to_string(),
                                sequence: results.template_sequence.clone(),
                            },
                            reference_file_name: String::new(),
                            reference_data: ReferenceData {
                                names: Vec::new(),
                                sequences: Vec::new(),
                            },
                            use_differential: results.differential_enabled,
                            exclusivity_file_names: Vec::new(),
                            exclusivity_data: None,
                            params: results.params.clone(),
                            output_folder: None,
                            template_length: results.template_length,
                            reference_count: results.total_sequences,
                            exclusivity_count: results
                                .exclusivity_sequence_count
                                .unwrap_or(0),
                        };
                        self.next_job_id += 1;

                        self.view_coverage_threshold = results.params.coverage_threshold;
                        self.differential_mode = results.differential_enabled;
                        self.results = Some(results.clone());
                        self.completed_jobs.push(CompletedJob { job, results });
                        self.selected_completed_job_index =
                            Some(self.completed_jobs.len() - 1);
                        self.load_error = None;
                        self.current_tab = Tab::Results;
                    }
                    Err(e) => {
                        self.load_error = Some(format!("Failed to parse: {}", e));
                    }
                },
                Err(e) => {
                    self.load_error = Some(format!("Failed to read file: {}", e));
                }
            }
        }
    }

    fn load_template_file(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("FASTA", &["fasta", "fa", "fna", "fas", "txt"])
            .pick_file()
        {
            match std::fs::read_to_string(&path) {
                Ok(content) => match parse_template_fasta(&content) {
                    Ok(data) => {
                        self.template_file_name = Some(
                            path.file_name()
                                .map(|n| n.to_string_lossy().to_string())
                                .unwrap_or_else(|| "unknown".to_string()),
                        );
                        self.template_data = Some(data);
                        self.template_error = None;
                    }
                    Err(e) => {
                        self.template_error = Some(e);
                    }
                },
                Err(e) => {
                    self.template_error = Some(format!("Failed to read file: {}", e));
                }
            }
        }
    }

    fn load_reference_file(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("FASTA", &["fasta", "fa", "fna", "fas", "txt"])
            .pick_file()
        {
            match std::fs::read_to_string(&path) {
                Ok(content) => match parse_reference_fasta(&content) {
                    Ok(data) => {
                        self.reference_file_name = Some(
                            path.file_name()
                                .map(|n| n.to_string_lossy().to_string())
                                .unwrap_or_else(|| "unknown".to_string()),
                        );
                        self.reference_data = Some(data);
                        self.reference_error = None;
                    }
                    Err(e) => {
                        self.reference_error = Some(e);
                    }
                },
                Err(e) => {
                    self.reference_error = Some(format!("Failed to read file: {}", e));
                }
            }
        }
    }

    fn add_exclusivity_file(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("FASTA", &["fasta", "fa", "fna", "fas", "txt"])
            .pick_file()
        {
            match std::fs::read_to_string(&path) {
                Ok(content) => match parse_reference_fasta(&content) {
                    Ok(data) => {
                        let file_name = path
                            .file_name()
                            .map(|n| n.to_string_lossy().to_string())
                            .unwrap_or_else(|| "unknown".to_string());
                        let min_len =
                            data.sequences.iter().map(|s| s.len()).min().unwrap_or(0);
                        let max_len =
                            data.sequences.iter().map(|s| s.len()).max().unwrap_or(0);
                        self.exclusivity_files.push(ExclusivityFileEntry {
                            file_name,
                            file_content: content,
                            sequence_count: data.len(),
                            min_length: min_len,
                            max_length: max_len,
                        });
                        self.rebuild_exclusivity_data();
                        self.exclusivity_error = None;
                    }
                    Err(e) => {
                        self.exclusivity_error = Some(e);
                    }
                },
                Err(e) => {
                    self.exclusivity_error = Some(format!("Failed to read file: {}", e));
                }
            }
        }
    }

    fn remove_exclusivity_file(&mut self, index: usize) {
        if index < self.exclusivity_files.len() {
            self.exclusivity_files.remove(index);
            self.rebuild_exclusivity_data();
        }
    }

    fn rebuild_exclusivity_data(&mut self) {
        if self.exclusivity_files.is_empty() {
            self.exclusivity_data = None;
            return;
        }

        let mut combined = ReferenceData::new();
        for entry in &self.exclusivity_files {
            if let Ok(data) = parse_reference_fasta(&entry.file_content) {
                combined.names.extend(data.names);
                combined.sequences.extend(data.sequences);
            }
        }

        if combined.sequences.is_empty() {
            self.exclusivity_data = None;
        } else {
            self.exclusivity_data = Some(combined);
        }
    }
}

impl AnalysisMethod {
    fn get_fixed_ambiguities(&self) -> u32 {
        match self {
            AnalysisMethod::FixedAmbiguities(n) => *n,
            _ => 1,
        }
    }

    fn get_incremental_pct(&self) -> u32 {
        match self {
            AnalysisMethod::Incremental(pct, _) => *pct,
            _ => 50,
        }
    }

    fn get_incremental_max_amb(&self) -> Option<u32> {
        match self {
            AnalysisMethod::Incremental(_, max_amb) => *max_amb,
            _ => None,
        }
    }
}

impl eframe::App for OligoscreenApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if self.is_analyzing {
            self.check_analysis_progress();
            ctx.request_repaint();
        }

        if self.pending_save {
            self.pending_save = false;
            self.save_results();
        }

        // Handle deferred exclusivity file removal
        if let Some(idx) = self.pending_remove_excl.take() {
            self.remove_exclusivity_file(idx);
        }

        // Top menu bar
        egui::TopBottomPanel::top("menu_bar").show(ctx, |ui| {
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Load Template...").clicked() {
                        self.load_template_file();
                        ui.close_menu();
                    }
                    if ui.button("Load References...").clicked() {
                        self.load_reference_file();
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui.button("Load Results from File...").clicked() {
                        self.load_results_into_completed();
                        ui.close_menu();
                    }
                    let can_save = self.results.is_some();
                    if ui
                        .add_enabled(can_save, egui::Button::new("Save Results..."))
                        .clicked()
                    {
                        self.save_results();
                        ui.close_menu();
                    }
                });
            });
        });

        // Tab bar
        egui::TopBottomPanel::top("tabs").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.current_tab, Tab::Input, "Input Data");
                ui.selectable_value(&mut self.current_tab, Tab::Analysis, "Analysis Setup");
                ui.selectable_value(
                    &mut self.current_tab,
                    Tab::Worklist,
                    format!("Worklist ({})", self.worklist.len()),
                );
                ui.selectable_value(
                    &mut self.current_tab,
                    Tab::Results,
                    format!("Results ({})", self.completed_jobs.len()),
                );
            });
        });

        // Status bar
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if self.is_analyzing {
                    ui.spinner();
                    if let Some(ref progress) = self.analysis_progress {
                        ui.label(format!("Processing: {}", &progress.message));
                    } else {
                        ui.label("Starting job...");
                    }
                } else if self.worklist_state == WorklistState::StopRequested {
                    ui.label("Stopping after current job...");
                } else {
                    let mut parts = Vec::new();
                    if !self.completed_jobs.is_empty() {
                        parts.push(format!(
                            "{} completed",
                            self.completed_jobs.len()
                        ));
                    }
                    if !self.worklist.is_empty() {
                        parts.push(format!("{} queued", self.worklist.len()));
                    }
                    if let Some(ref t) = self.template_data {
                        parts.push(format!("Template: {} bp", t.sequence.len()));
                    }
                    if let Some(ref r) = self.reference_data {
                        parts.push(format!("References: {} seqs", r.len()));
                    }
                    if parts.is_empty() {
                        ui.label("Load template and reference sequences to begin");
                    } else {
                        ui.label(parts.join(" | "));
                    }
                }
            });
        });

        // Main content
        egui::CentralPanel::default().show(ctx, |ui| {
            match self.current_tab {
                Tab::Input => self.show_input_tab(ui),
                Tab::Analysis => self.show_analysis_tab(ui),
                Tab::Worklist => self.show_worklist_tab(ui),
                Tab::Results => self.show_results_tab(ui),
            }
        });

        // Detail window
        if self.show_detail_window {
            self.show_variant_detail_window(ctx);
        }
    }
}

impl OligoscreenApp {
    fn show_input_tab(&mut self, ui: &mut egui::Ui) {
        ui.heading("Input Data");
        ui.separator();

        // --- Template Sequence ---
        ui.group(|ui| {
            ui.horizontal(|ui| {
                ui.heading("Template Sequence");
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.button("Clear").clicked() {
                        self.template_file_name = None;
                        self.template_data = None;
                        self.template_error = None;
                    }
                    if ui.button("Load File").clicked() {
                        self.load_template_file();
                    }
                });
            });

            ui.label("Single sequence in FASTA format (A, C, G, T only)");

            if let Some(ref error) = self.template_error {
                ui.colored_label(egui::Color32::RED, format!("Error: {}", error));
            }
            if let Some(ref data) = self.template_data {
                ui.horizontal(|ui| {
                    ui.colored_label(
                        egui::Color32::from_rgb(100, 200, 100),
                        format!(
                            "File: {}",
                            self.template_file_name.as_deref().unwrap_or("unknown")
                        ),
                    );
                });
                ui.colored_label(
                    egui::Color32::from_rgb(100, 200, 100),
                    format!("Sequence: {} ({} bp)", data.name, data.sequence.len()),
                );
            } else {
                ui.colored_label(egui::Color32::GRAY, "No template loaded");
            }
        });

        ui.add_space(5.0);

        // --- Reference Sequences ---
        ui.group(|ui| {
            ui.horizontal(|ui| {
                ui.heading("Reference Sequences");
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.button("Clear").clicked() {
                        self.reference_file_name = None;
                        self.reference_data = None;
                        self.reference_error = None;
                    }
                    if ui.button("Load File").clicked() {
                        self.load_reference_file();
                    }
                });
            });

            ui.label("Multiple sequences in FASTA format (unaligned)");

            if let Some(ref error) = self.reference_error {
                ui.colored_label(egui::Color32::RED, format!("Error: {}", error));
            }
            if let Some(ref data) = self.reference_data {
                let min_len = data.sequences.iter().map(|s| s.len()).min().unwrap_or(0);
                let max_len = data.sequences.iter().map(|s| s.len()).max().unwrap_or(0);
                ui.horizontal(|ui| {
                    ui.colored_label(
                        egui::Color32::from_rgb(100, 200, 100),
                        format!(
                            "File: {}",
                            self.reference_file_name.as_deref().unwrap_or("unknown")
                        ),
                    );
                });
                ui.colored_label(
                    egui::Color32::from_rgb(100, 200, 100),
                    format!(
                        "{} sequences ({}-{} bp)",
                        data.len(),
                        min_len,
                        max_len
                    ),
                );
            } else {
                ui.colored_label(egui::Color32::GRAY, "No references loaded");
            }
        });

        ui.add_space(10.0);

        // --- Differential Analysis / Exclusivity Sequences ---
        ui.checkbox(&mut self.use_differential, "Use differential analysis");

        if self.use_differential {
            ui.group(|ui| {
                ui.horizontal(|ui| {
                    ui.heading("Exclusivity Sequences");
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.button("Add File").clicked() {
                            self.add_exclusivity_file();
                        }
                    });
                });

                ui.label("Sequences that oligos must be distinct from (off-targets)");

                if let Some(ref error) = self.exclusivity_error {
                    ui.colored_label(egui::Color32::RED, format!("Error: {}", error));
                }

                if self.exclusivity_files.is_empty() {
                    ui.colored_label(egui::Color32::GRAY, "No exclusivity files loaded");
                } else {
                    let mut remove_idx = None;
                    for (i, entry) in self.exclusivity_files.iter().enumerate() {
                        ui.horizontal(|ui| {
                            if ui.small_button("X").clicked() {
                                remove_idx = Some(i);
                            }
                            ui.label(format!(
                                "{} - {} sequences ({}-{} bp)",
                                entry.file_name,
                                entry.sequence_count,
                                entry.min_length,
                                entry.max_length
                            ));
                        });
                    }
                    if let Some(idx) = remove_idx {
                        self.pending_remove_excl = Some(idx);
                    }

                    // Summary
                    if let Some(ref data) = self.exclusivity_data {
                        ui.separator();
                        ui.colored_label(
                            egui::Color32::from_rgb(100, 200, 100),
                            format!(
                                "Total: {} exclusivity sequences from {} file(s)",
                                data.len(),
                                self.exclusivity_files.len()
                            ),
                        );
                    }
                }
            });
        }

        ui.add_space(10.0);

        // --- Output Folder ---
        ui.group(|ui| {
            ui.horizontal(|ui| {
                ui.heading("Output Folder (Optional)");
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if self.output_folder.is_some() {
                        if ui.button("Clear").clicked() {
                            self.output_folder = None;
                        }
                    }
                    if ui.button("Select Folder").clicked() {
                        self.select_output_folder();
                    }
                });
            });
            ui.label(
                "If set, results will be auto-saved as JSON to this folder after analysis.",
            );
            if let Some(ref folder) = self.output_folder {
                ui.colored_label(egui::Color32::from_rgb(100, 200, 100), format!("Folder: {}", folder));
            } else {
                ui.colored_label(egui::Color32::GRAY, "No output folder selected (manual save only)");
            }
        });

        ui.add_space(10.0);

        // --- Add to Worklist ---
        let can_add = self.template_data.is_some() && self.reference_data.is_some();
        let warn_excl =
            self.use_differential && self.exclusivity_data.is_none();
        ui.horizontal(|ui| {
            if ui
                .add_enabled(can_add, egui::Button::new("Add to Worklist"))
                .clicked()
            {
                self.add_to_worklist();
            }
            if !can_add {
                ui.colored_label(
                    egui::Color32::GRAY,
                    "Load template and references first",
                );
            }
            if warn_excl {
                ui.colored_label(
                    egui::Color32::YELLOW,
                    "Differential enabled but no exclusivity files loaded",
                );
            }
        });
    }

    fn show_analysis_tab(&mut self, ui: &mut egui::Ui) {
        ui.heading("Analysis Setup");
        ui.separator();
        ui.label("These settings apply to all jobs added to the worklist.");
        ui.add_space(5.0);

        egui::ScrollArea::vertical().show(ui, |ui| {
            // Pairwise Aligner Settings
            ui.group(|ui| {
                ui.heading("Pairwise Aligner Settings");

                ui.horizontal(|ui| {
                    ui.label("Match score:");
                    ui.add(
                        egui::DragValue::new(&mut self.params.pairwise.match_score).range(0..=10),
                    );
                    ui.add_space(20.0);
                    ui.label("Mismatch score:");
                    ui.add(
                        egui::DragValue::new(&mut self.params.pairwise.mismatch_score)
                            .range(-10..=0),
                    );
                });

                ui.horizontal(|ui| {
                    ui.label("Gap open penalty:");
                    ui.add(
                        egui::DragValue::new(&mut self.params.pairwise.gap_open_penalty)
                            .range(-20..=0),
                    );
                    ui.add_space(20.0);
                    ui.label("Gap extend penalty:");
                    ui.add(
                        egui::DragValue::new(&mut self.params.pairwise.gap_extend_penalty)
                            .range(-20..=0),
                    );
                });

                ui.horizontal(|ui| {
                    ui.label("Maximum allowed mismatches:");
                    ui.add(
                        egui::DragValue::new(&mut self.params.pairwise.max_mismatches)
                            .range(0..=50),
                    );
                });
                ui.label("Matches exceeding this mismatch count are recorded as 'no match'.");
            });

            ui.add_space(10.0);

            // Analysis method selection
            ui.group(|ui| {
                ui.heading("Analysis Method");

                ui.radio_value(
                    &mut self.method_selection,
                    MethodSelection::NoAmbiguities,
                    "No Ambiguities - Find all unique exact variants",
                );

                ui.horizontal(|ui| {
                    ui.radio_value(
                        &mut self.method_selection,
                        MethodSelection::FixedAmbiguities,
                        "Fixed Ambiguities - Use up to N ambiguity codes per variant",
                    );
                });

                if self.method_selection == MethodSelection::FixedAmbiguities {
                    ui.horizontal(|ui| {
                        ui.add_space(20.0);
                        ui.label("Max ambiguities:");
                        let mut n = self.params.method.get_fixed_ambiguities();
                        if ui.add(egui::DragValue::new(&mut n).range(0..=20)).changed() {
                            self.params.method = AnalysisMethod::FixedAmbiguities(n);
                        }
                    });
                }

                ui.horizontal(|ui| {
                    ui.radio_value(
                        &mut self.method_selection,
                        MethodSelection::Incremental,
                        "Incremental - Find variants covering X% of remaining sequences",
                    );
                });

                if self.method_selection == MethodSelection::Incremental {
                    ui.horizontal(|ui| {
                        ui.add_space(20.0);
                        ui.label("Target coverage per step (%):");
                        let mut pct = self.params.method.get_incremental_pct();
                        let max_amb = self.params.method.get_incremental_max_amb();
                        if ui
                            .add(egui::DragValue::new(&mut pct).range(1..=100))
                            .changed()
                        {
                            self.params.method = AnalysisMethod::Incremental(pct, max_amb);
                        }
                    });
                    ui.horizontal(|ui| {
                        ui.add_space(20.0);
                        ui.checkbox(
                            &mut self.incremental_limit_ambiguities,
                            "Limit ambiguities:",
                        );
                        ui.add_enabled(
                            self.incremental_limit_ambiguities,
                            egui::DragValue::new(&mut self.incremental_max_ambiguities)
                                .range(0..=20),
                        );
                        ui.label("max");
                    });
                    if self.incremental_limit_ambiguities {
                        ui.horizontal(|ui| {
                            ui.add_space(20.0);
                            ui.label(
                                "If target % cannot be reached, accepts best variant within limit.",
                            );
                        });
                    }
                }
            });

            ui.add_space(10.0);

            // Global options
            ui.group(|ui| {
                ui.heading("Global Options");
                ui.checkbox(
                    &mut self.params.exclude_n,
                    "Exclude N (any base) as ambiguity code",
                );
            });

            ui.add_space(10.0);

            // Oligo length range
            ui.group(|ui| {
                ui.heading("Oligo Length Range");
                ui.horizontal(|ui| {
                    ui.label("Minimum length:");
                    ui.add(
                        egui::DragValue::new(&mut self.params.min_oligo_length).range(3..=100),
                    );
                    ui.add_space(20.0);
                    ui.label("Maximum length:");
                    ui.add(
                        egui::DragValue::new(&mut self.params.max_oligo_length).range(3..=100),
                    );
                });

                if self.params.min_oligo_length > self.params.max_oligo_length {
                    self.params.max_oligo_length = self.params.min_oligo_length;
                }

                let range = self.params.max_oligo_length - self.params.min_oligo_length + 1;
                if range > 20 {
                    ui.colored_label(
                        egui::Color32::YELLOW,
                        format!(
                            "Warning: Large length range ({}) may take significant time",
                            range
                        ),
                    );
                }
            });

            ui.add_space(10.0);

            // Resolution
            ui.group(|ui| {
                ui.heading("Analysis Resolution");
                ui.horizontal(|ui| {
                    ui.label("Step size (bases):");
                    ui.add(egui::DragValue::new(&mut self.params.resolution).range(1..=100));
                });
                ui.label("Lower values = more positions analyzed, higher resolution");
            });

            ui.add_space(10.0);

            // Coverage threshold
            ui.group(|ui| {
                ui.heading("Coverage Threshold");
                ui.horizontal(|ui| {
                    ui.label("Target coverage (%):");
                    ui.add(
                        egui::DragValue::new(&mut self.params.coverage_threshold)
                            .range(1.0..=100.0),
                    );
                });
                ui.label("Number of variants needed to reach this coverage will be reported");
            });

        });
    }

    fn show_worklist_tab(&mut self, ui: &mut egui::Ui) {
        ui.heading("Worklist");
        ui.separator();

        // === Parallelization (moved from Analysis Setup) ===
        ui.group(|ui| {
            ui.heading("Parallelization");

            let available_threads = std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1);

            ui.label(format!("Available parallelism: {} threads", available_threads));

            ui.horizontal(|ui| {
                ui.radio_value(
                    &mut self.thread_selection,
                    ThreadSelection::Auto,
                    format!("Auto ({} threads)", available_threads),
                );
            });
            ui.horizontal(|ui| {
                ui.radio_value(
                    &mut self.thread_selection,
                    ThreadSelection::Manual,
                    "Manual:",
                );
                let enabled = self.thread_selection == ThreadSelection::Manual;
                ui.add_enabled(
                    enabled,
                    egui::DragValue::new(&mut self.manual_thread_count)
                        .range(1..=available_threads.max(32)),
                );
                ui.label("threads");
            });
        });

        ui.add_space(10.0);

        // === Process / Stop Controls ===
        ui.horizontal(|ui| {
            let can_process =
                !self.worklist.is_empty() && self.worklist_state == WorklistState::Idle;
            if ui
                .add_enabled(can_process, egui::Button::new("Process Worklist"))
                .clicked()
            {
                self.start_worklist_processing();
            }

            let can_stop = self.worklist_state == WorklistState::Processing;
            if ui
                .add_enabled(can_stop, egui::Button::new("Stop After Current"))
                .clicked()
            {
                self.worklist_state = WorklistState::StopRequested;
            }

            match self.worklist_state {
                WorklistState::Idle => {}
                WorklistState::Processing => {
                    ui.spinner();
                    let jobs_done =
                        self.worklist_total_at_start - self.worklist.len();
                    ui.label(format!(
                        "Processing job {} of {}",
                        jobs_done + 1,
                        self.worklist_total_at_start
                    ));
                }
                WorklistState::StopRequested => {
                    ui.spinner();
                    ui.colored_label(
                        egui::Color32::YELLOW,
                        "Stopping after current job...",
                    );
                }
            }
        });

        ui.add_space(5.0);

        // === Progress Bars ===
        if self.worklist_state != WorklistState::Idle {
            let jobs_done = self.worklist_total_at_start - self.worklist.len();
            let overall_frac = if self.worklist_total_at_start > 0 {
                jobs_done as f32 / self.worklist_total_at_start as f32
            } else {
                0.0
            };
            ui.horizontal(|ui| {
                ui.label("Overall:");
                ui.add(
                    egui::ProgressBar::new(overall_frac).text(format!(
                        "{}/{} jobs",
                        jobs_done, self.worklist_total_at_start
                    )),
                );
            });

            if let Some(ref progress) = self.analysis_progress {
                let job_frac = if progress.total_lengths > 0 {
                    let length_frac =
                        progress.lengths_completed as f32 / progress.total_lengths as f32;
                    let pos_frac = if progress.total_positions > 0 {
                        // Use completed count from the message (parsed from "Position X/Y")
                        // Fall back to a rough estimate from position index
                        (progress.lengths_completed as f32
                            + (1.0 / progress.total_lengths as f32))
                            .min(1.0)
                    } else {
                        0.0
                    };
                    let _ = pos_frac;
                    length_frac
                } else {
                    0.0
                };
                ui.horizontal(|ui| {
                    ui.label("Current job:");
                    ui.add(
                        egui::ProgressBar::new(job_frac).text(&progress.message),
                    );
                });
            }
        }

        ui.add_space(10.0);

        // === Queued Jobs Table ===
        ui.heading("Queued Jobs");
        if self.worklist.is_empty() {
            ui.colored_label(
                egui::Color32::GRAY,
                "No jobs queued. Use the Input Data tab to add jobs.",
            );
        } else {
            let mut pending_remove: Option<usize> = None;

            egui::ScrollArea::vertical()
                .id_salt("worklist_scroll")
                .max_height(300.0)
                .show(ui, |ui| {
                    egui::Grid::new("worklist_grid")
                        .striped(true)
                        .min_col_width(40.0)
                        .show(ui, |ui| {
                            // Header
                            ui.strong("");
                            ui.strong("#");
                            ui.strong("Template");
                            ui.strong("References");
                            ui.strong("Exclusivity");
                            ui.strong("Oligo Range");
                            ui.strong("Method");
                            ui.strong("Output");
                            ui.end_row();

                            for (i, job) in self.worklist.iter().enumerate() {
                                let is_current =
                                    self.worklist_state == WorklistState::Processing
                                        && i == self.current_job_index;

                                if is_current {
                                    ui.spinner();
                                } else if ui.small_button("X").clicked() {
                                    pending_remove = Some(i);
                                }

                                ui.label(format!("{}", job.id));
                                ui.label(&job.template_file_name);
                                ui.label(format!("{} seqs", job.reference_count));
                                if job.use_differential {
                                    ui.label(format!("{} seqs", job.exclusivity_count));
                                } else {
                                    ui.label("-");
                                }
                                ui.label(format!(
                                    "{}-{} bp",
                                    job.params.min_oligo_length,
                                    job.params.max_oligo_length
                                ));
                                ui.label(job.params.method.description());
                                if job.output_folder.is_some() {
                                    ui.label("Auto-save");
                                } else {
                                    ui.label("-");
                                }
                                ui.end_row();
                            }
                        });
                });

            if let Some(idx) = pending_remove {
                self.remove_worklist_job(idx);
            }
        }

        // === Completed Jobs Summary ===
        if !self.completed_jobs.is_empty() {
            ui.add_space(10.0);
            ui.separator();
            ui.label(format!(
                "{} completed job(s) available in the Results tab.",
                self.completed_jobs.len()
            ));
        }

        // === Auto-save error ===
        if let Some(ref err) = self.auto_save_error {
            ui.colored_label(egui::Color32::RED, err);
        }
    }

    fn show_results_tab(&mut self, ui: &mut egui::Ui) {
        if self.completed_jobs.is_empty() {
            ui.heading("Results");
            ui.separator();
            ui.label(
                "No completed jobs yet. Add jobs in the Input tab and process them in the Worklist tab.",
            );
            ui.add_space(10.0);
            if ui.button("Load Results from File").clicked() {
                self.load_results_into_completed();
            }
            if let Some(ref error) = self.load_error {
                ui.colored_label(egui::Color32::RED, error);
            }
            return;
        }

        // Job selector + header
        ui.horizontal(|ui| {
            ui.heading("Results");

            ui.separator();
            ui.label("Job:");

            let selected_label = self
                .selected_completed_job_index
                .and_then(|i| self.completed_jobs.get(i))
                .map(|cj| {
                    format!("#{} - {}", cj.job.id, cj.job.template_file_name)
                })
                .unwrap_or_else(|| "Select a job".to_string());

            let mut new_selection = self.selected_completed_job_index;
            egui::ComboBox::from_id_salt("completed_job_selector")
                .selected_text(&selected_label)
                .show_ui(ui, |ui| {
                    for (i, cj) in self.completed_jobs.iter().enumerate() {
                        let label = format!(
                            "#{} - {} ({} refs, {}-{} bp)",
                            cj.job.id,
                            cj.job.template_file_name,
                            cj.job.reference_count,
                            cj.job.params.min_oligo_length,
                            cj.job.params.max_oligo_length,
                        );
                        ui.selectable_value(&mut new_selection, Some(i), label);
                    }
                });

            // Sync results when selection changes
            if new_selection != self.selected_completed_job_index {
                self.selected_completed_job_index = new_selection;
                if let Some(idx) = new_selection {
                    if let Some(cj) = self.completed_jobs.get(idx) {
                        self.results = Some(cj.results.clone());
                        self.view_coverage_threshold = cj.results.params.coverage_threshold;
                        self.differential_mode = cj.results.differential_enabled;
                    }
                }
            }

            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui.button("Load Results from File").clicked() {
                    self.load_results_into_completed();
                }
                let has_results = self.results.is_some();
                if ui
                    .add_enabled(has_results, egui::Button::new("Save Results"))
                    .clicked()
                {
                    self.pending_save = true;
                }
            });
        });
        ui.separator();

        if self.results.is_none() {
            ui.label("Select a completed job to view its results.");
            return;
        }

        // Extract data we need
        let (lengths, template_seq, total_seqs, has_differential) = {
            let results = self.results.as_ref().unwrap();
            let mut lengths: Vec<u32> = results.results_by_length.keys().copied().collect();
            lengths.sort();
            (
                lengths,
                results.template_sequence.clone(),
                results.total_sequences,
                results.differential_enabled,
            )
        };

        if lengths.is_empty() {
            ui.label("No length results available.");
            return;
        }

        // Controls row 1: zoom + info + differential toggle
        ui.horizontal(|ui| {
            ui.label("Zoom:");
            ui.add(egui::Slider::new(&mut self.zoom_level, 0.5..=3.0));
            ui.add_space(20.0);
            ui.label(format!(
                "{} reference sequences | Template: {} bp",
                total_seqs,
                template_seq.len()
            ));
            if has_differential {
                ui.separator();
                ui.checkbox(&mut self.differential_mode, "Differential mode");
            }
        });

        if !self.differential_mode {
            // === NORMAL MODE CONTROLS ===

            // Controls row 2: coverage threshold + color range
            ui.horizontal(|ui| {
                ui.label("Coverage threshold (%):");
                ui.add(
                    egui::DragValue::new(&mut self.view_coverage_threshold)
                        .range(1.0..=100.0)
                        .speed(0.5),
                );
                if ui.button("Apply").clicked() {
                    self.recalculate_coverage_threshold();
                }
                ui.separator();
                ui.label("Color range - Green at:");
                ui.add(egui::DragValue::new(&mut self.color_green_at).range(1..=1000));
                ui.label("variants, Red at:");
                ui.add(egui::DragValue::new(&mut self.color_red_at).range(1..=1000));
                ui.label("variants");
            });

            // Ensure green <= red
            if self.color_green_at > self.color_red_at {
                self.color_red_at = self.color_green_at;
            }

            // Controls row 3: no-match darkening thresholds
            ui.horizontal(|ui| {
                ui.label("No-match darkening - OK at:");
                ui.add(
                    egui::DragValue::new(&mut self.nomatch_ok_percent)
                        .range(0.0..=100.0)
                        .speed(0.5)
                        .suffix("%"),
                );
                ui.label(", Dark red at:");
                ui.add(
                    egui::DragValue::new(&mut self.nomatch_bad_percent)
                        .range(0.0..=100.0)
                        .speed(0.5)
                        .suffix("%"),
                );
            });

            if self.nomatch_ok_percent > self.nomatch_bad_percent {
                self.nomatch_bad_percent = self.nomatch_ok_percent;
            }
        } else {
            // === DIFFERENTIAL MODE CONTROLS ===

            // Exclusivity color controls
            ui.horizontal(|ui| {
                ui.label("Exclusivity color - Green at:");
                ui.add(egui::DragValue::new(&mut self.diff_green_at).range(0..=50));
                ui.label("mismatches, Red at:");
                ui.add(egui::DragValue::new(&mut self.diff_red_at).range(0..=50));
                ui.label("mismatches");
                ui.separator();
                ui.label("Ignore best:");
                ui.add(egui::DragValue::new(&mut self.diff_ignore_count).range(0..=1000));
                ui.label("sequences");
            });

            // Darkening controls (conservation metrics)
            ui.horizontal(|ui| {
                ui.label("Darkening - Variant count: Green at:");
                ui.add(egui::DragValue::new(&mut self.color_green_at).range(1..=1000));
                ui.label(", Red at:");
                ui.add(egui::DragValue::new(&mut self.color_red_at).range(1..=1000));
                ui.separator();
                ui.label("No-match: OK at:");
                ui.add(
                    egui::DragValue::new(&mut self.nomatch_ok_percent)
                        .range(0.0..=100.0)
                        .speed(0.5)
                        .suffix("%"),
                );
                ui.label(", Bad at:");
                ui.add(
                    egui::DragValue::new(&mut self.nomatch_bad_percent)
                        .range(0.0..=100.0)
                        .speed(0.5)
                        .suffix("%"),
                );
            });

            if self.color_green_at > self.color_red_at {
                self.color_red_at = self.color_green_at;
            }
            if self.nomatch_ok_percent > self.nomatch_bad_percent {
                self.nomatch_bad_percent = self.nomatch_ok_percent;
            }

            // Coverage threshold (still needed for variant count)
            ui.horizontal(|ui| {
                ui.label("Coverage threshold (%):");
                ui.add(
                    egui::DragValue::new(&mut self.view_coverage_threshold)
                        .range(1.0..=100.0)
                        .speed(0.5),
                );
                if ui.button("Apply").clicked() {
                    self.recalculate_coverage_threshold();
                }
            });
        }

        ui.add_space(5.0);

        // Heatmap display
        let coverage_threshold = self.view_coverage_threshold;
        self.show_heatmap(ui, &lengths, &template_seq, coverage_threshold);

        // Error messages
        if let Some(ref error) = self.save_error {
            ui.colored_label(egui::Color32::RED, error);
        }
        if let Some(ref error) = self.load_error {
            ui.colored_label(egui::Color32::RED, error);
        }
    }

    fn show_heatmap(
        &mut self,
        ui: &mut egui::Ui,
        lengths: &[u32],
        template_seq: &str,
        coverage_threshold: f64,
    ) {
        let results = self.results.as_ref().unwrap();

        // Get positions from the first length result
        let first_length_result = results.results_by_length.get(&lengths[0]);
        let positions: Vec<usize> = first_length_result
            .map(|lr| lr.positions.iter().map(|p| p.position).collect())
            .unwrap_or_default();

        if positions.is_empty() {
            ui.label("No positions analyzed.");
            return;
        }

        // Cell dimensions: zoom only affects horizontal width, height is fixed
        let cell_w = (14.0 * self.zoom_level).max(3.0);
        let cell_h: f32 = 54.0;
        let label_width: f32 = 50.0;
        let header_height: f32 = 20.0;
        let pos_label_height: f32 = 14.0;

        let num_cols = positions.len();
        let num_rows = lengths.len();

        // Summary stats per length
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                for &length in lengths {
                    if let Some(lr) = results.results_by_length.get(&length) {
                        let non_skipped: Vec<_> =
                            lr.positions.iter().filter(|p| !p.analysis.skipped).collect();
                        if !non_skipped.is_empty() {
                            let avg: f64 =
                                non_skipped.iter().map(|p| p.variants_needed).sum::<usize>()
                                    as f64
                                    / non_skipped.len() as f64;
                            let min = non_skipped
                                .iter()
                                .map(|p| p.variants_needed)
                                .min()
                                .unwrap_or(0);
                            let max = non_skipped
                                .iter()
                                .map(|p| p.variants_needed)
                                .max()
                                .unwrap_or(0);
                            ui.label(format!(
                                "{}bp: {}-{} (avg {:.1})",
                                length, min, max, avg
                            ));
                            ui.separator();
                        }
                    }
                }
            });
        });

        ui.add_space(5.0);

        if self.differential_mode {
            ui.label("Exclusivity: min mismatches (green=specific, red=similar to off-targets). Darkened by conservation metrics.");
        } else {
            ui.label(format!(
                "Variants needed to reach {:.0}% coverage (click cell for details):",
                coverage_threshold
            ));
        }

        // Build heatmap data: lookup by (length, position)
        let heatmap_data: std::collections::HashMap<
            (u32, usize),
            &crate::analysis::PositionResult,
        > = {
            let mut map = std::collections::HashMap::new();
            for &length in lengths {
                if let Some(lr) = results.results_by_length.get(&length) {
                    for pr in &lr.positions {
                        map.insert((length, pr.position), pr);
                    }
                }
            }
            map
        };

        // Total width/height for the heatmap area
        let total_width = label_width + (num_cols as f32 * cell_w);
        let total_height =
            pos_label_height + header_height + (num_rows as f32 * cell_h) + 30.0;

        let scroll_output = egui::ScrollArea::horizontal()
            .id_salt("heatmap_scroll")
            .show(ui, |ui| {
                let (response, painter) = ui.allocate_painter(
                    egui::vec2(total_width, total_height),
                    egui::Sense::click_and_drag(),
                );
                let origin = response.rect.min;

                // --- Position numbers row ---
                let show_every_n = if cell_w < 12.0 {
                    (12.0 / cell_w).ceil() as usize
                } else {
                    1
                };

                for (col, &pos) in positions.iter().enumerate() {
                    if col % show_every_n != 0 {
                        continue;
                    }
                    let x = origin.x + label_width + (col as f32 * cell_w) + cell_w / 2.0;
                    let y = origin.y + pos_label_height / 2.0;
                    painter.text(
                        egui::pos2(x, y),
                        egui::Align2::CENTER_CENTER,
                        format!("{}", pos + 1),
                        egui::FontId::proportional(9.0),
                        egui::Color32::GRAY,
                    );
                }

                // --- Template sequence row ---
                let seq_y_start = origin.y + pos_label_height;
                if cell_w >= 8.0 {
                    for (col, &pos) in positions.iter().enumerate() {
                        if pos < template_seq.len() {
                            let base = &template_seq[pos..pos + 1];
                            let x =
                                origin.x + label_width + (col as f32 * cell_w) + cell_w / 2.0;
                            let y = seq_y_start + header_height / 2.0;

                            let color = base_color(base.chars().next().unwrap_or('N'));
                            painter.text(
                                egui::pos2(x, y),
                                egui::Align2::CENTER_CENTER,
                                base,
                                egui::FontId::monospace(11.0),
                                color,
                            );
                        }
                    }
                } else {
                    for (col, &pos) in positions.iter().enumerate() {
                        if pos < template_seq.len() {
                            let base_char = template_seq.as_bytes()[pos] as char;
                            let color = base_color(base_char);
                            let x = origin.x + label_width + (col as f32 * cell_w);
                            let tick_rect = egui::Rect::from_min_size(
                                egui::pos2(x, seq_y_start + 2.0),
                                egui::vec2((cell_w - 1.0).max(1.0), header_height - 4.0),
                            );
                            painter.rect_filled(tick_rect, 0.0, color);
                        }
                    }
                }

                // --- Row labels (oligo lengths) ---
                let grid_y_start = seq_y_start + header_height;
                for (row, &length) in lengths.iter().enumerate() {
                    let y = grid_y_start + (row as f32 * cell_h) + cell_h / 2.0;
                    painter.text(
                        egui::pos2(origin.x + label_width - 5.0, y),
                        egui::Align2::RIGHT_CENTER,
                        format!("{} bp", length),
                        egui::FontId::proportional(11.0),
                        egui::Color32::LIGHT_GRAY,
                    );
                }

                // --- Heatmap cells ---
                let mut hovered_cell: Option<(u32, usize)> = None;
                let mut clicked_cell: Option<(u32, usize)> = None;

                let is_differential = self.differential_mode;

                for (row, &length) in lengths.iter().enumerate() {
                    for (col, &pos) in positions.iter().enumerate() {
                        let cell_x = origin.x + label_width + (col as f32 * cell_w);
                        let cell_y = grid_y_start + (row as f32 * cell_h);
                        let cell_rect = egui::Rect::from_min_size(
                            egui::pos2(cell_x, cell_y),
                            egui::vec2(cell_w - 1.0, cell_h - 1.0),
                        );

                        let color = if let Some(pr) = heatmap_data.get(&(length, pos)) {
                            if pr.analysis.skipped {
                                egui::Color32::from_rgb(40, 40, 40)
                            } else if is_differential {
                                let eff_min_mm = pr
                                    .exclusivity
                                    .as_ref()
                                    .map(|e| {
                                        effective_min_mismatches(e, self.diff_ignore_count)
                                    })
                                    .flatten();
                                let no_match_frac = if pr.analysis.total_sequences > 0 {
                                    pr.analysis.no_match_count as f64
                                        / pr.analysis.total_sequences as f64
                                } else {
                                    0.0
                                };
                                differential_position_color(
                                    eff_min_mm,
                                    pr.variants_needed,
                                    no_match_frac,
                                    self.diff_green_at,
                                    self.diff_red_at,
                                    self.color_green_at,
                                    self.color_red_at,
                                    self.nomatch_ok_percent / 100.0,
                                    self.nomatch_bad_percent / 100.0,
                                )
                            } else {
                                let no_match_frac = if pr.analysis.total_sequences > 0 {
                                    pr.analysis.no_match_count as f64
                                        / pr.analysis.total_sequences as f64
                                } else {
                                    0.0
                                };
                                position_color(
                                    pr.variants_needed,
                                    no_match_frac,
                                    self.color_green_at,
                                    self.color_red_at,
                                    self.nomatch_ok_percent / 100.0,
                                    self.nomatch_bad_percent / 100.0,
                                )
                            }
                        } else {
                            egui::Color32::from_rgb(30, 30, 30)
                        };

                        painter.rect_filled(cell_rect, 1.0, color);

                        if let Some(pointer_pos) = response.hover_pos() {
                            if cell_rect.contains(pointer_pos) {
                                hovered_cell = Some((length, pos));
                                painter.rect_stroke(
                                    cell_rect,
                                    1.0,
                                    egui::Stroke::new(1.5, egui::Color32::WHITE),
                                    egui::StrokeKind::Outside,
                                );
                            }
                        }

                        if response.clicked() {
                            if let Some(pointer_pos) = ui.ctx().pointer_latest_pos() {
                                if cell_rect.contains(pointer_pos) {
                                    clicked_cell = Some((length, pos));
                                }
                            }
                        }
                    }
                }

                // Handle tooltip
                if let Some((length, pos)) = hovered_cell {
                    if let Some(pr) = heatmap_data.get(&(length, pos)) {
                        let mut tooltip_text = if pr.analysis.skipped {
                            format!(
                                "Position: {}, Length: {} bp\nSkipped: {}",
                                pos + 1,
                                length,
                                pr.analysis
                                    .skip_reason
                                    .as_deref()
                                    .unwrap_or("Unknown")
                            )
                        } else {
                            format!(
                                "Position: {}, Length: {} bp\nVariants needed: {}\nCoverage: {:.1}%\nMatched: {}/{}\nNo match: {}",
                                pos + 1,
                                length,
                                pr.variants_needed,
                                pr.analysis.coverage_at_threshold,
                                pr.analysis.sequences_analyzed,
                                pr.analysis.total_sequences,
                                pr.analysis.no_match_count,
                            )
                        };

                        // Add exclusivity info to tooltip
                        if let Some(ref excl) = pr.exclusivity {
                            let eff = effective_min_mismatches(excl, self.diff_ignore_count);
                            let mm_str = match eff {
                                Some(mm) => format!("{}", mm),
                                None => "all no-match".to_string(),
                            };
                            tooltip_text.push_str(&format!(
                                "\nExclusivity: min mismatches = {} ({} sequences)",
                                mm_str, excl.total_sequences
                            ));
                        }

                        response.clone().on_hover_text(tooltip_text);
                    }
                }

                // Handle click
                if let Some((length, pos)) = clicked_cell {
                    self.selected_position = Some(pos);
                    self.selected_length_for_detail = Some(length);
                    self.show_detail_window = true;
                }
            });

        // Redirect vertical mouse wheel to horizontal scroll when hovering over heatmap
        if let Some(hover_pos) = ui.ctx().pointer_hover_pos() {
            if scroll_output.inner_rect.contains(hover_pos) {
                let vertical_delta = ui.input(|i| i.smooth_scroll_delta.y);
                if vertical_delta.abs() > 0.1 {
                    let mut state = scroll_output.state;
                    state.offset.x -= vertical_delta;
                    state.offset.x = state.offset.x.clamp(
                        0.0,
                        (total_width - scroll_output.inner_rect.width()).max(0.0),
                    );
                    state.store(ui.ctx(), scroll_output.id);
                    ui.ctx().request_repaint();
                }
            }
        }

        // Legend
        ui.add_space(5.0);
        if self.differential_mode {
            self.show_differential_legend(ui);
        } else {
            self.show_normal_legend(ui);
        }
    }

    fn show_normal_legend(&self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label("Legend:");
            ui.add_space(10.0);

            let g = self.color_green_at;
            let r = self.color_red_at;
            let sample_points: Vec<(usize, String)> = if r <= g {
                vec![(g, format!("<={}", g)), (g + 1, format!(">{}", g))]
            } else {
                let mid = (g + r) / 2;
                let mut pts = vec![(g, format!("<={}", g))];
                if mid > g && mid < r {
                    pts.push((mid, format!("{}", mid)));
                }
                pts.push((r, format!(">={}", r)));
                pts
            };

            let nm_ok = self.nomatch_ok_percent / 100.0;
            let nm_bad = self.nomatch_bad_percent / 100.0;

            for (count, label) in &sample_points {
                let color = position_color(*count, 0.0, g, r, nm_ok, nm_bad);
                let (rect, _) =
                    ui.allocate_exact_size(egui::vec2(15.0, 15.0), egui::Sense::hover());
                ui.painter().rect_filled(rect, 2.0, color);
                ui.label(label);
                ui.add_space(8.0);
            }

            ui.separator();

            let mid_count = (g + r) / 2;
            let mid_count = if mid_count < 1 { 1 } else { mid_count };
            let nm_samples = [
                (nm_ok, format!("{}%", self.nomatch_ok_percent as u32)),
                (nm_bad, format!("{}%", self.nomatch_bad_percent as u32)),
            ];
            ui.label("No-match:");
            for (nm_frac, label) in &nm_samples {
                let color = position_color(mid_count, *nm_frac, g, r, nm_ok, nm_bad);
                let (rect, _) =
                    ui.allocate_exact_size(egui::vec2(15.0, 15.0), egui::Sense::hover());
                ui.painter().rect_filled(rect, 2.0, color);
                ui.label(label);
                ui.add_space(4.0);
            }

            ui.separator();
            let (rect, _) =
                ui.allocate_exact_size(egui::vec2(15.0, 15.0), egui::Sense::hover());
            ui.painter()
                .rect_filled(rect, 2.0, egui::Color32::from_rgb(40, 40, 40));
            ui.label("skipped/no data");
        });
    }

    fn show_differential_legend(&self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label("Legend (Differential):");
            ui.add_space(10.0);

            // Exclusivity gradient samples (no darkening)
            let dg = self.diff_green_at;
            let dr = self.diff_red_at;

            let sample_mms: Vec<(Option<u32>, String)> = if dg > dr {
                vec![
                    (Some(dg), format!(">={} mm", dg)),
                    (Some((dg + dr) / 2), format!("{} mm", (dg + dr) / 2)),
                    (Some(dr), format!("<={} mm", dr)),
                ]
            } else {
                vec![
                    (Some(dg), format!("{} mm", dg)),
                    (Some(dr), format!("{} mm", dr)),
                ]
            };

            for (mm_val, label) in &sample_mms {
                let color = differential_position_color(
                    *mm_val, 1, 0.0, dg, dr, self.color_green_at, self.color_red_at, 1.0, 1.0,
                );
                let (rect, _) =
                    ui.allocate_exact_size(egui::vec2(15.0, 15.0), egui::Sense::hover());
                ui.painter().rect_filled(rect, 2.0, color);
                ui.label(label);
                ui.add_space(4.0);
            }

            ui.separator();
            ui.label("+ darkening from conservation");

            ui.separator();
            let (rect, _) =
                ui.allocate_exact_size(egui::vec2(15.0, 15.0), egui::Sense::hover());
            ui.painter()
                .rect_filled(rect, 2.0, egui::Color32::from_rgb(40, 40, 40));
            ui.label("skipped/no data");
        });
    }

    fn show_variant_detail_window(&mut self, ctx: &egui::Context) {
        let Some(ref results) = self.results else {
            self.show_detail_window = false;
            return;
        };

        let Some(length) = self.selected_length_for_detail else {
            self.show_detail_window = false;
            return;
        };

        let Some(position) = self.selected_position else {
            self.show_detail_window = false;
            return;
        };

        let Some(length_result) = results.results_by_length.get(&length) else {
            self.show_detail_window = false;
            return;
        };

        let Some(pos_result) = length_result
            .positions
            .iter()
            .find(|p| p.position == position)
        else {
            self.show_detail_window = false;
            return;
        };

        let pos_result = pos_result.clone();
        let coverage_threshold = results.params.coverage_threshold;

        // Extract template oligo for display
        let template_oligo = if position + length as usize <= results.template_sequence.len() {
            &results.template_sequence[position..position + length as usize]
        } else {
            ""
        };
        let template_oligo = template_oligo.to_string();

        let show_reverse_complement = self.detail_show_reverse_complement;
        let show_codon_spacing = self.detail_show_codon_spacing;

        egui::Window::new(format!("Position {} Details", position + 1))
            .open(&mut self.show_detail_window)
            .default_width(650.0)
            .default_height(500.0)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.label(format!("Position: {}", position + 1));
                    ui.separator();
                    ui.label(format!("Oligo length: {} bp", length));
                });

                // Template oligo display
                if !template_oligo.is_empty() {
                    let display_template = format_sequence_for_display(
                        &template_oligo,
                        show_reverse_complement,
                        show_codon_spacing,
                    );
                    ui.horizontal(|ui| {
                        ui.label("Template oligo:");
                        ui.add(
                            egui::Label::new(
                                egui::RichText::new(&display_template)
                                    .monospace()
                                    .size(11.0)
                                    .color(egui::Color32::from_rgb(100, 180, 255)),
                            )
                            .wrap_mode(egui::TextWrapMode::Extend),
                        );
                    });
                }

                ui.separator();

                if pos_result.analysis.skipped {
                    ui.colored_label(
                        egui::Color32::YELLOW,
                        format!(
                            "This window was skipped: {}",
                            pos_result
                                .analysis
                                .skip_reason
                                .as_deref()
                                .unwrap_or("Unknown reason")
                        ),
                    );
                    return;
                }

                ui.label(format!(
                    "Total references: {}",
                    pos_result.analysis.total_sequences
                ));
                ui.label(format!(
                    "Matched: {}",
                    pos_result.analysis.sequences_analyzed
                ));
                if pos_result.analysis.no_match_count > 0 {
                    ui.colored_label(
                        egui::Color32::from_rgb(255, 180, 100),
                        format!(
                            "No match: {}/{} ({:.1}%)",
                            pos_result.analysis.no_match_count,
                            pos_result.analysis.total_sequences,
                            (pos_result.analysis.no_match_count as f64
                                / pos_result.analysis.total_sequences as f64)
                                * 100.0
                        ),
                    );
                }
                ui.label(format!(
                    "Variants needed for {:.0}% coverage: {}",
                    coverage_threshold, pos_result.variants_needed
                ));
                ui.label(format!(
                    "Coverage at threshold: {:.1}%",
                    pos_result.analysis.coverage_at_threshold
                ));

                ui.separator();

                // Display options
                ui.horizontal(|ui| {
                    ui.heading("Variants");
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        ui.checkbox(&mut self.detail_show_codon_spacing, "Codon spacing");
                        ui.checkbox(
                            &mut self.detail_show_reverse_complement,
                            "Reverse complement",
                        );
                    });
                });

                egui::ScrollArea::vertical()
                    .id_salt("detail_scroll")
                    .max_height(250.0)
                    .show(ui, |ui| {
                        egui::Grid::new("variants_grid")
                            .striped(true)
                            .min_col_width(50.0)
                            .show(ui, |ui| {
                                ui.strong("#");
                                ui.strong("Sequence");
                                ui.strong("Count");
                                ui.strong("Percentage");
                                ui.strong("Cumulative");
                                ui.end_row();

                                let mut cumulative = 0.0;
                                for (i, variant) in
                                    pos_result.analysis.variants.iter().enumerate()
                                {
                                    cumulative += variant.percentage;

                                    let is_threshold = i + 1 == pos_result.variants_needed;

                                    if is_threshold {
                                        ui.colored_label(
                                            egui::Color32::GREEN,
                                            format!("{}", i + 1),
                                        );
                                    } else {
                                        ui.label(format!("{}", i + 1));
                                    }

                                    let display_seq = format_sequence_for_display(
                                        &variant.sequence,
                                        show_reverse_complement,
                                        show_codon_spacing,
                                    );

                                    ui.add(
                                        egui::Label::new(
                                            egui::RichText::new(&display_seq)
                                                .monospace()
                                                .size(11.0),
                                        )
                                        .wrap_mode(egui::TextWrapMode::Extend),
                                    );

                                    ui.label(format!("{}", variant.count));
                                    ui.label(format!("{:.1}%", variant.percentage));

                                    if is_threshold {
                                        ui.colored_label(
                                            egui::Color32::GREEN,
                                            format!("{:.1}%", cumulative),
                                        );
                                    } else {
                                        ui.label(format!("{:.1}%", cumulative));
                                    }

                                    ui.end_row();
                                }

                                // No match row
                                if pos_result.analysis.no_match_count > 0 {
                                    ui.label("");
                                    ui.colored_label(
                                        egui::Color32::from_rgb(255, 180, 100),
                                        "No match",
                                    );
                                    ui.colored_label(
                                        egui::Color32::from_rgb(255, 180, 100),
                                        format!("{}", pos_result.analysis.no_match_count),
                                    );
                                    let no_match_pct = (pos_result.analysis.no_match_count
                                        as f64
                                        / pos_result.analysis.total_sequences as f64)
                                        * 100.0;
                                    ui.colored_label(
                                        egui::Color32::from_rgb(255, 180, 100),
                                        format!("{:.1}%", no_match_pct),
                                    );
                                    ui.label("");
                                    ui.end_row();
                                }
                            });

                        // === Exclusivity Analysis Section ===
                        if let Some(ref excl) = pos_result.exclusivity {
                            ui.add_space(10.0);
                            ui.separator();
                            ui.heading("Exclusivity Analysis");

                            ui.label(format!(
                                "Total exclusivity sequences: {}",
                                excl.total_sequences
                            ));
                            if let Some(min_mm) = excl.min_mismatches {
                                ui.label(format!("Minimum mismatches: {}", min_mm));
                            } else {
                                ui.colored_label(
                                    egui::Color32::from_rgb(100, 200, 100),
                                    "All exclusivity sequences: no match (fully specific)",
                                );
                            }

                            ui.add_space(5.0);

                            egui::Grid::new("exclusivity_grid")
                                .striped(true)
                                .min_col_width(60.0)
                                .show(ui, |ui| {
                                    ui.strong("Mismatches");
                                    ui.strong("Count");
                                    ui.strong("Example");
                                    ui.end_row();

                                    for bucket in &excl.mismatch_histogram {
                                        if bucket.mismatches == u32::MAX {
                                            ui.colored_label(
                                                egui::Color32::from_rgb(100, 200, 100),
                                                "No match",
                                            );
                                        } else {
                                            let color = if bucket.mismatches == 0 {
                                                egui::Color32::from_rgb(255, 80, 80)
                                            } else if bucket.mismatches <= 2 {
                                                egui::Color32::from_rgb(255, 180, 100)
                                            } else {
                                                egui::Color32::LIGHT_GRAY
                                            };
                                            ui.colored_label(
                                                color,
                                                format!("{}", bucket.mismatches),
                                            );
                                        }
                                        ui.label(format!("{}", bucket.count));
                                        ui.label(&bucket.example_name);
                                        ui.end_row();
                                    }
                                });
                        }
                    });
            });
    }
}

/// Calculate effective minimum mismatches after ignoring the best N sequences.
fn effective_min_mismatches(
    excl: &crate::analysis::ExclusivityResult,
    ignore_count: usize,
) -> Option<u32> {
    if ignore_count == 0 {
        return excl.min_mismatches;
    }

    let mut remaining_ignore = ignore_count;
    for bucket in &excl.mismatch_histogram {
        if bucket.mismatches == u32::MAX {
            // No-match bucket — these are already "infinite", skip them
            continue;
        }
        if bucket.count <= remaining_ignore {
            remaining_ignore -= bucket.count;
        } else {
            // This bucket has sequences remaining after ignoring
            return Some(bucket.mismatches);
        }
    }

    // All matched sequences were ignored — effectively all are no-match
    None
}

/// Format a sequence for display with optional transformations
fn format_sequence_for_display(seq: &str, reverse_comp: bool, codon_spacing: bool) -> String {
    let mut result = if reverse_comp {
        reverse_complement(seq)
    } else {
        seq.to_string()
    };

    if codon_spacing {
        result = add_codon_spacing(&result);
    }

    result
}

/// Add spaces every 3 characters (codon format)
fn add_codon_spacing(seq: &str) -> String {
    seq.chars()
        .enumerate()
        .flat_map(|(i, c)| {
            if i > 0 && i % 3 == 0 {
                vec![' ', c]
            } else {
                vec![c]
            }
        })
        .collect()
}

/// Get color for a position based on variant count and no-match fraction (normal mode).
fn position_color(
    variant_count: usize,
    no_match_fraction: f64,
    green_at: usize,
    red_at: usize,
    nomatch_ok: f64,
    nomatch_bad: f64,
) -> egui::Color32 {
    if variant_count == 0 {
        return egui::Color32::from_rgb(40, 40, 40);
    }

    let (base_r, base_g, base_b) =
        green_yellow_red_gradient(variant_count, green_at, red_at);

    // No-match darkening
    let dark_red = (100.0f64, 20.0f64, 20.0f64);
    let nm_t = ramp(no_match_fraction, nomatch_ok, nomatch_bad);

    let r = (base_r * (1.0 - nm_t) + dark_red.0 * nm_t).clamp(0.0, 255.0) as u8;
    let g = (base_g * (1.0 - nm_t) + dark_red.1 * nm_t).clamp(0.0, 255.0) as u8;
    let b = (base_b * (1.0 - nm_t) + dark_red.2 * nm_t).clamp(0.0, 255.0) as u8;

    egui::Color32::from_rgb(r, g, b)
}

/// Get color for a position in differential mode.
///
/// Base color: exclusivity min mismatches gradient (green=high=specific, red=low=similar).
/// Darkening: conservation metrics (variant count + no-match %) blend toward dark red.
fn differential_position_color(
    min_mismatches: Option<u32>,
    variant_count: usize,
    no_match_fraction: f64,
    diff_green_at: u32,
    diff_red_at: u32,
    var_green_at: usize,
    var_red_at: usize,
    nomatch_ok: f64,
    nomatch_bad: f64,
) -> egui::Color32 {
    // Conservation darkening always applies — compute it first.
    // If either metric reaches its worst threshold, the cell goes fully dark red
    // regardless of how good the exclusivity score is.
    let variant_dark = ramp_usize(variant_count, var_green_at, var_red_at);
    let nomatch_dark = ramp(no_match_fraction, nomatch_ok, nomatch_bad);
    let darkening = variant_dark.max(nomatch_dark);

    // Skipped positions (zero variants analyzed) → dark gray
    if variant_count == 0 {
        return egui::Color32::from_rgb(40, 40, 40);
    }

    // Base color from exclusivity: green→yellow→red gradient
    // None = all no-match = fully specific = best = green (t=0)
    let t = match min_mismatches {
        None => 0.0,
        Some(mm) => {
            if diff_green_at <= diff_red_at {
                if mm <= diff_green_at { 0.0 } else { 1.0 }
            } else if mm >= diff_green_at {
                0.0
            } else if mm <= diff_red_at {
                1.0
            } else {
                (diff_green_at - mm) as f64 / (diff_green_at - diff_red_at) as f64
            }
        }
    };

    let (base_r, base_g, base_b) = green_yellow_red_from_t(t);

    // Blend base color toward dark red by the darkening factor
    let dark_red = (100.0f64, 20.0f64, 20.0f64);
    let r = (base_r * (1.0 - darkening) + dark_red.0 * darkening).clamp(0.0, 255.0) as u8;
    let g = (base_g * (1.0 - darkening) + dark_red.1 * darkening).clamp(0.0, 255.0) as u8;
    let b = (base_b * (1.0 - darkening) + dark_red.2 * darkening).clamp(0.0, 255.0) as u8;

    egui::Color32::from_rgb(r, g, b)
}

/// 3-stop gradient: green → yellow → red. Returns (r, g, b) as f64.
fn green_yellow_red_gradient(value: usize, green_at: usize, red_at: usize) -> (f64, f64, f64) {
    let t = if red_at <= green_at {
        if value <= green_at {
            0.0
        } else {
            1.0
        }
    } else if value <= green_at {
        0.0
    } else if value >= red_at {
        1.0
    } else {
        (value - green_at) as f64 / (red_at - green_at) as f64
    };

    green_yellow_red_from_t(t)
}

/// Convert t (0..1) to green→yellow→red gradient RGB.
fn green_yellow_red_from_t(t: f64) -> (f64, f64, f64) {
    let green = (0.0f64, 180.0f64, 0.0f64);
    let yellow = (220.0f64, 200.0f64, 0.0f64);
    let red = (220.0f64, 50.0f64, 50.0f64);

    if t <= 0.5 {
        let s = t * 2.0;
        (
            green.0 + (yellow.0 - green.0) * s,
            green.1 + (yellow.1 - green.1) * s,
            green.2 + (yellow.2 - green.2) * s,
        )
    } else {
        let s = (t - 0.5) * 2.0;
        (
            yellow.0 + (red.0 - yellow.0) * s,
            yellow.1 + (red.1 - yellow.1) * s,
            yellow.2 + (red.2 - yellow.2) * s,
        )
    }
}

/// Convert t=0 to green color (for "all no-match" case in differential mode).
fn green_yellow_red_to_color(t: f64) -> egui::Color32 {
    let (r, g, b) = green_yellow_red_from_t(t);
    egui::Color32::from_rgb(r as u8, g as u8, b as u8)
}

/// Linear ramp: 0 at low, 1 at high, clamped.
fn ramp(value: f64, low: f64, high: f64) -> f64 {
    let v = value.clamp(0.0, 1.0);
    let lo = low.clamp(0.0, 1.0);
    let hi = high.clamp(0.0, 1.0);
    if hi <= lo {
        if v <= lo {
            0.0
        } else {
            1.0
        }
    } else if v <= lo {
        0.0
    } else if v >= hi {
        1.0
    } else {
        (v - lo) / (hi - lo)
    }
}

/// Linear ramp for usize values.
fn ramp_usize(value: usize, low: usize, high: usize) -> f64 {
    if high <= low {
        if value <= low {
            0.0
        } else {
            1.0
        }
    } else if value <= low {
        0.0
    } else if value >= high {
        1.0
    } else {
        (value - low) as f64 / (high - low) as f64
    }
}

/// Color for DNA base letters in the template display
fn base_color(base: char) -> egui::Color32 {
    match base {
        'A' => egui::Color32::from_rgb(100, 200, 100), // Green
        'T' => egui::Color32::from_rgb(220, 80, 80),   // Red
        'G' => egui::Color32::from_rgb(255, 200, 60),   // Yellow/gold
        'C' => egui::Color32::from_rgb(100, 150, 255),  // Blue
        _ => egui::Color32::GRAY,
    }
}
