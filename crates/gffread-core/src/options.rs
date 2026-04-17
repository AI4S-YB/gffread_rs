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
    pub cds: Option<PathBuf>,
    pub protein: Option<PathBuf>,
    pub write_exon_segments: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RuntimeOptions {
    pub expose_warnings: bool,
    pub output: Option<PathBuf>,
    pub main_output: MainOutput,
    pub table_format: Option<String>,
    pub genome: Option<PathBuf>,
    pub fasta_outputs: FastaOutputs,
    pub input: PathBuf,
    pub inputs: Vec<PathBuf>,
    pub original_args: Vec<String>,
}
