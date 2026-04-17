# GffRead C++ to Rust Function Map

This document records the phase-one compatibility surface implemented in this worktree. The local C++ `gffread` binary remains the oracle for observable behavior; the Rust modules below are the current ports or compatibility shims for that surface.

## Process and CLI

- `gffread.cpp::main` maps to [`crates/gffread-cli/src/main.rs`](../../crates/gffread-cli/src/main.rs), [`crates/gffread-cli/src/parse.rs`](../../crates/gffread-cli/src/parse.rs), and [`crates/gffread-cli/src/process.rs`](../../crates/gffread-cli/src/process.rs).
- Command-line mode dispatch (`--version`, help/error flow, run mode) maps to `parse::parse_args`, `parse::CommandMode`, and `process::run_process`.
- Upstream output-file opening paths, including the phase-one "main output plus FASTA outputs" wiring, map to `process::run_outputs`.
- `gffread.cpp::printGff3Header` maps to `gffread_core::emit::gff3::write_gff3`.
- Upstream command-line echoing in `-E` mode maps to `process::oracle_command_line` plus the warning block emitted from `process::run_outputs`.

## Parsing and Model

- The subset of GCLib `GffLoader::load` / `GffReader` behavior needed by the phase-one example fixtures maps to `gffread_core::loader::gff::load_annotation` and `parse_annotation`.
- The phase-one transcript record model maps from `GffObj` to `gffread_core::model::Transcript`.
- Exon and CDS segment storage maps from GCLib exon/CDS coordinate lists to `gffread_core::model::Segment`.
- Runtime option state that was previously held across `gffread.cpp` globals/locals maps to `gffread_core::options::RuntimeOptions`, `MainOutput`, and `FastaOutputs`.
- Compatibility-formatted process errors map to `gffread_core::compat::CompatError`.

## Output

- GFF3 transcript/exon/CDS emission needed by the current oracle suite maps to `gffread_core::emit::gff3::write_gff3` and its internal transcript writer.
- GFF3 locus feature emission for the current `-M`/`--cluster-only` oracle suite maps to `gffread_core::emit::locus::write_locus`.
- GTF transcript/exon/CDS emission needed by the current oracle suite maps to `gffread_core::emit::gtf::write_gtf`.
- `gffread.cpp::setTableFormat` maps to `gffread_core::emit::table::parse_format`.
- `gff_utils.cpp::printTableData` maps to `gffread_core::emit::table::write_table`.
- Upstream file/stdout routing for GFF3, GTF, and table output maps to the `MainOutput` branch handling in `crates/gffread-cli/src/process.rs`.

## Filtering and Clustering

- `loadIDlist` plus `GffLoader::checkFilters` ID include/exclude behavior maps to `gffread_core::filters::apply_filters` and the `IdFilter` runtime options.
- `GffLoader::checkFilters` range, minimum covered length, maximum intron, no-pseudo, and attribute filtering behavior maps to `gffread_core::filters::apply_filters`, `RangeFilter`, and the GFF3 attribute emission controls.
- The current oracle-covered subset of `gff_utils.cpp::collectLocusData` maps to `gffread_core::cluster::apply_clustering` and `gffread_core::model::Locus`.
- The current oracle-covered subset of `GffLoader::redundantTranscripts` maps to `gffread_core::cluster::collapse_redundancy`, including exact intron-chain redundancy, `-Q` relaxed boundary containment, and `-K` contained intron-chain collapse for the minimized fixtures.

## Sequence Logic

- The FASTA loading/index-writing behavior needed by the phase-one examples maps to `gffread_core::fasta::load_genome`.
- Upstream spliced transcript extraction maps to `gffread_core::fasta::spliced_sequence` and `write_transcript_fasta`.
- Phase-one CDS extraction with leading/trailing phase trimming maps to `gffread_core::fasta::spliced_cds_sequence`, `trimmed_segment_bounds`, and `write_cds_fasta`.
- FASTA defline/record writing maps to `gffread_core::fasta::write_fasta_record`, `transcript_defline`, and `projected_defline`.
- Reverse-complement logic maps to `gffread_core::fasta::reverse_complement`.
- Codon translation and protein FASTA emission map to `gffread_core::fasta::translate`, `translate_codon`, and `write_protein_fasta`.

## Phase-One Boundary

- `gff_utils.cpp::adjust_stopcodon` remains for a later compatibility phase.
- Full `collectLocusData` parity beyond the current minimized `-M`/`--cluster-only` fixtures remains for a later compatibility phase.
- Full `GffLoader::redundantTranscripts` parity beyond exact intron-chain, `-Q`, and `-K` minimized fixtures remains for a later compatibility phase.
- Advanced validation, splice/CDS filters, streaming, BED/TLF, complex locus metadata, and the broader non-phase-one CLI surface remain for later phases.
