use std::fs::File;
use std::io::{self, Write};
use std::path::PathBuf;

use gffread_core::compat::CompatError;
use gffread_core::emit::{gff3, gtf, table};
use gffread_core::loader::gff::load_annotation;
use gffread_core::options::{MainOutput, RuntimeOptions};
use gffread_core::VERSION;

use crate::help::USAGE;
use crate::parse::{parse_args, CommandMode};

pub fn run_process(args: Vec<String>) -> i32 {
    match parse_args(args) {
        Ok(CommandMode::Version) => {
            print!("{VERSION}\n");
            0
        }
        Ok(CommandMode::Help) => {
            let _ = io::stderr().write_all(USAGE.as_bytes());
            1
        }
        Ok(CommandMode::Run(options)) => match run_outputs(options) {
            Ok(()) => 0,
            Err(err) => {
                let _ = io::stderr().write_all(err.message.as_bytes());
                err.exit_code
            }
        },
        Err(err) => {
            let _ = io::stderr().write_all(err.message.as_bytes());
            err.exit_code
        }
    }
}

fn run_outputs(options: RuntimeOptions) -> Result<(), CompatError> {
    if options.genome.is_some()
        || options.fasta_outputs.transcript.is_some()
        || options.fasta_outputs.cds.is_some()
        || options.fasta_outputs.protein.is_some()
        || options.fasta_outputs.write_exon_segments
    {
        return Err(CompatError::new(
            "Error: FASTA-related options are not implemented in this phase-one step\n",
            1,
        ));
    }

    let annotation = load_annotation(&options.input)?;
    let command_line = oracle_command_line(&options.original_args);

    if options.expose_warnings {
        eprint!(
            "Command line was:\n{command_line}\n   .. loaded {} genomic features from {}\n",
            annotation.transcripts.len(),
            options.input.display()
        );
    }

    match options.main_output {
        MainOutput::Gff3 => {
            let output = options.output.ok_or_else(|| {
                CompatError::new("Error: output file is required in phase-one GFF3 path\n", 1)
            })?;
            let mut file = File::create(&output).map_err(|_| {
                CompatError::new(format!("Error creating file: {}\n", output.display()), 1)
            })?;

            gff3::write_gff3(&mut file, &annotation, VERSION, &command_line)
                .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
            Ok(())
        }
        MainOutput::Gtf => {
            let output = options.output.ok_or_else(|| {
                CompatError::new("Error: output file is required in phase-one GTF path\n", 1)
            })?;
            let mut file = File::create(&output).map_err(|_| {
                CompatError::new(format!("Error creating file: {}\n", output.display()), 1)
            })?;

            gtf::write_gtf(&mut file, &annotation)
                .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
            Ok(())
        }
        MainOutput::Table => {
            let output = options.output.ok_or_else(|| {
                CompatError::new(
                    "Error: output file is required in phase-one table path\n",
                    1,
                )
            })?;
            let format = options
                .table_format
                .as_deref()
                .ok_or_else(|| CompatError::new("Error: --table requires a format\n", 1))?;
            let mut file = File::create(&output).map_err(|_| {
                CompatError::new(format!("Error creating file: {}\n", output.display()), 1)
            })?;

            table::write_table(&mut file, &annotation, format)
                .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
            Ok(())
        }
    }
}

fn oracle_command_line(args: &[String]) -> String {
    let oracle = workspace_root().join("gffread");
    if args.is_empty() {
        oracle.display().to_string()
    } else {
        format!("{} {}", oracle.display(), args.join(" "))
    }
}

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .canonicalize()
        .expect("workspace root must resolve")
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::run_outputs;
    use gffread_core::options::{FastaOutputs, MainOutput, RuntimeOptions};

    fn example_path(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("..")
            .join("..")
            .join("examples")
            .join(name)
    }

    fn temp_output_path(case_name: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock should be after Unix epoch")
            .as_nanos();
        std::env::temp_dir().join(format!(
            "gffread-rs-{case_name}-{}-{nanos}.gff",
            std::process::id()
        ))
    }

    fn gff3_options(
        output: &Path,
        fasta_outputs: FastaOutputs,
        genome: Option<PathBuf>,
    ) -> RuntimeOptions {
        RuntimeOptions {
            expose_warnings: false,
            output: Some(output.to_path_buf()),
            main_output: MainOutput::Gff3,
            table_format: None,
            genome,
            fasta_outputs,
            input: example_path("annotation.gff"),
            inputs: vec![example_path("annotation.gff")],
            original_args: Vec::new(),
        }
    }

    #[test]
    fn fasta_related_options_are_rejected_before_gff3_output() {
        let cases = [
            (
                "genome",
                Some(example_path("genome.fa")),
                FastaOutputs::default(),
            ),
            (
                "transcript_fasta",
                None,
                FastaOutputs {
                    transcript: Some(PathBuf::from("transcripts.fa")),
                    ..FastaOutputs::default()
                },
            ),
            (
                "cds_fasta",
                None,
                FastaOutputs {
                    cds: Some(PathBuf::from("transcripts_CDS.fa")),
                    ..FastaOutputs::default()
                },
            ),
            (
                "protein_fasta",
                None,
                FastaOutputs {
                    protein: Some(PathBuf::from("transcripts_prot.fa")),
                    ..FastaOutputs::default()
                },
            ),
            (
                "exon_segments",
                None,
                FastaOutputs {
                    write_exon_segments: true,
                    ..FastaOutputs::default()
                },
            ),
        ];

        for (case_name, genome, fasta_outputs) in cases {
            let output = temp_output_path(case_name);
            let result = run_outputs(gff3_options(&output, fasta_outputs, genome));

            let err =
                result.expect_err("FASTA-related options must not silently write GFF3 output");
            assert_eq!(err.exit_code, 1);
            assert!(
                err.message.contains("FASTA-related options"),
                "unexpected error message: {}",
                err.message
            );
            assert!(
                !output.exists(),
                "GFF3 output should not be created for {case_name}"
            );

            let _ = fs::remove_file(output);
        }
    }
}
