use std::path::PathBuf;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MainOutput {
    Gff3,
    Gtf,
    Table,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FastaOutputs {
    pub transcript: Option<PathBuf>,
    pub unspliced: Option<PathBuf>,
    pub cds: Option<PathBuf>,
    pub protein: Option<PathBuf>,
    pub write_exon_segments: bool,
    pub suppress_transcript_cds: bool,
    pub write_protein_star_stop: bool,
    pub padding: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IdFilterMode {
    Include,
    Exclude,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RangeFilter {
    pub seqid: String,
    pub strand: Option<char>,
    pub start: u64,
    pub end: u64,
    pub fully_within: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IdFilter {
    pub path: PathBuf,
    pub mode: IdFilterMode,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ClusterOptions {
    pub merge: bool,
    pub cluster_only: bool,
    pub relax_boundary_containment: bool,
    pub collapse_contained: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RuntimeOptions {
    pub expose_warnings: bool,
    pub output: Option<PathBuf>,
    pub main_output: MainOutput,
    pub table_format: Option<String>,
    pub attrs: Option<Vec<String>>,
    pub keep_all_attrs: bool,
    pub genome: Option<PathBuf>,
    pub fasta_outputs: FastaOutputs,
    pub range_filter: Option<RangeFilter>,
    pub id_filter: Option<IdFilter>,
    pub min_length: Option<u64>,
    pub max_intron: Option<u64>,
    pub coding_only: bool,
    pub noncoding_only: bool,
    pub multi_exon_only: bool,
    pub no_pseudo: bool,
    pub cluster: ClusterOptions,
    pub input: PathBuf,
    pub inputs: Vec<PathBuf>,
    pub original_args: Vec<String>,
}
