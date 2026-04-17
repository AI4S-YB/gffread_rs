use std::fs::File;
use std::io::{self, Write};
use std::path::PathBuf;

use gffread_core::compat::CompatError;
use gffread_core::emit::{gff3, gtf, table};
use gffread_core::fasta::{load_genome, write_cds_fasta, write_transcript_fasta};
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
    if options.fasta_outputs.protein.is_some() {
        return Err(CompatError::new(
            "Error: protein FASTA output (-y) is not implemented in this phase-one step\n",
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
            if let Some(output) = &options.output {
                let mut file = File::create(output).map_err(|_| {
                    CompatError::new(format!("Error creating file: {}\n", output.display()), 1)
                })?;

                gff3::write_gff3(&mut file, &annotation, VERSION, &command_line)
                    .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
            }
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
        }
    }

    if options.fasta_outputs.transcript.is_some() || options.fasta_outputs.cds.is_some() {
        let genome_path = options.genome.as_ref().ok_or_else(|| {
            CompatError::new("Error: -g option is required for options -w/x/y/u/V/N/M !\n", 1)
        })?;
        let genome = load_genome(genome_path)?;

        if let Some(path) = &options.fasta_outputs.transcript {
            let mut file = File::create(path).map_err(|_| {
                CompatError::new(format!("Error creating file: {}\n", path.display()), 1)
            })?;
            write_transcript_fasta(&mut file, &annotation, &genome)?;
        }

        if let Some(path) = &options.fasta_outputs.cds {
            let mut file = File::create(path).map_err(|_| {
                CompatError::new(format!("Error creating file: {}\n", path.display()), 1)
            })?;
            write_cds_fasta(
                &mut file,
                &annotation,
                &genome,
                options.fasta_outputs.write_exon_segments,
            )?;
        }
    }

    Ok(())
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

    struct TempFixtures {
        root: PathBuf,
        annotation: PathBuf,
        genome: PathBuf,
    }

    impl Drop for TempFixtures {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.root);
        }
    }

    fn copied_example_fixtures(case_name: &str) -> TempFixtures {
        let root = std::env::temp_dir().join(format!(
            "gffread-rs-fixtures-{case_name}-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("system clock should be after Unix epoch")
                .as_nanos()
        ));
        fs::create_dir_all(&root).expect("temporary fixture directory should be created");

        let annotation = root.join("annotation.gff");
        fs::copy(example_path("annotation.gff"), &annotation)
            .expect("annotation fixture should be copied");

        let genome = root.join("genome.fa");
        fs::copy(example_path("genome.fa"), &genome).expect("genome fixture should be copied");

        TempFixtures {
            root,
            annotation,
            genome,
        }
    }

    fn gff3_options(
        input: &Path,
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
            input: input.to_path_buf(),
            inputs: vec![input.to_path_buf()],
            original_args: Vec::new(),
        }
    }

    #[test]
    fn protein_fasta_option_is_still_rejected() {
        let fixtures = copied_example_fixtures("protein_fasta");
        let output = temp_output_path("protein_fasta");
        let result = run_outputs(gff3_options(
            &fixtures.annotation,
            &output,
            FastaOutputs {
                protein: Some(PathBuf::from("transcripts_prot.fa")),
                ..FastaOutputs::default()
            },
            Some(fixtures.genome.clone()),
        ));

        let err = result.expect_err("protein FASTA should remain unimplemented");
        assert_eq!(err.exit_code, 1);
        assert!(
            err.message.contains("protein FASTA output"),
            "unexpected error message: {}",
            err.message
        );
        assert!(!output.exists(), "GFF3 output should not be created");

        let _ = fs::remove_file(output);
    }

    #[test]
    fn transcript_and_cds_fasta_options_can_run_with_gff3_output() {
        let fixtures = copied_example_fixtures("transcript_cds_fasta");
        let output = temp_output_path("transcript_cds_fasta");
        let transcript_fasta = std::env::temp_dir().join(format!(
            "gffread-rs-transcripts-{}-{}.fa",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("system clock should be after Unix epoch")
                .as_nanos()
        ));
        let cds_fasta = std::env::temp_dir().join(format!(
            "gffread-rs-cds-{}-{}.fa",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("system clock should be after Unix epoch")
                .as_nanos()
        ));
        let fai_path = example_path("genome.fa.fai");

        let result = run_outputs(gff3_options(
            &fixtures.annotation,
            &output,
            FastaOutputs {
                transcript: Some(transcript_fasta.clone()),
                cds: Some(cds_fasta.clone()),
                write_exon_segments: true,
                ..FastaOutputs::default()
            },
            Some(fixtures.genome.clone()),
        ));

        result.expect("transcript and CDS FASTA should be supported");
        assert!(output.exists(), "GFF3 output should still be written");
        assert!(transcript_fasta.exists(), "transcript FASTA should be written");
        assert!(cds_fasta.exists(), "CDS FASTA should be written");
        assert!(
            !fai_path.exists(),
            "process tests must not create FASTA indices in tracked examples"
        );

        let _ = fs::remove_file(output);
        let _ = fs::remove_file(transcript_fasta);
        let _ = fs::remove_file(cds_fasta);
    }
}
