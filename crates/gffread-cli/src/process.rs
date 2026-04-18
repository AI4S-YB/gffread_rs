use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

use gffread_core::cluster::apply_clustering;
use gffread_core::compat::CompatError;
use gffread_core::emit::{bed, tlf};
use gffread_core::emit::{gff3, gtf, table};
use gffread_core::fasta::{
    load_genome, write_cds_fasta, write_protein_fasta, write_transcript_fasta,
    write_unspliced_fasta,
};
use gffread_core::filters::apply_filters;
use gffread_core::loader::gff::load_annotation;
use gffread_core::model::Annotation;
use gffread_core::options::{MainOutput, RuntimeOptions};
use gffread_core::sort::{sort_annotation, sort_transcripts_for_fasta};
use gffread_core::VERSION;

use crate::help::USAGE;
use crate::parse::{parse_args, CommandMode};

pub fn run_process(program: String, args: Vec<String>) -> i32 {
    match parse_args(program, args) {
        Ok(CommandMode::Version) => {
            println!("{VERSION}");
            0
        }
        Ok(CommandMode::Help) => {
            let _ = io::stderr().write_all(USAGE.as_bytes());
            1
        }
        Ok(CommandMode::UsageError(err)) => {
            if err.show_usage {
                let _ = io::stderr().write_all(USAGE.as_bytes());
                let _ = io::stderr().write_all(b"\n");
            }
            let _ = io::stderr().write_all(err.message.as_bytes());
            err.exit_code
        }
        Ok(CommandMode::Run(options)) => match run_outputs(*options) {
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
    let need_genome = options.fasta_outputs.transcript.is_some()
        || options.fasta_outputs.unspliced.is_some()
        || options.fasta_outputs.cds.is_some()
        || options.fasta_outputs.protein.is_some();
    let genome = if need_genome {
        let genome_path = options.genome.as_ref().ok_or_else(|| {
            CompatError::new(
                "Error: -g option is required for options -w/x/y/u/V/N/M !\n",
                1,
            )
        })?;
        Some(load_genome(genome_path)?)
    } else {
        None
    };
    let annotation = load_annotation(&options.input, options.input_format)?;
    let filtered_annotation = apply_filters(&annotation, &options)?;
    let mut annotation = filtered_annotation;
    let command_line = command_line(&options.program, &options.original_args);

    if !options.cluster.merge {
        let fasta_annotation = Annotation {
            transcripts: sort_transcripts_for_fasta(
                &annotation.transcripts,
                &options.ref_sort_order,
                &annotation.ref_order,
            )?,
            genes: Vec::new(),
            ref_order: annotation.ref_order.clone(),
            header_comments: Vec::new(),
        };
        write_fasta_outputs(&options, &fasta_annotation, genome.as_ref())?;
    }

    sort_annotation(&mut annotation, &options.ref_sort_order, options.keep_genes)?;
    let (mut annotation, mut loci) = apply_clustering(&annotation, &options.cluster);
    if options.cluster.merge {
        write_fasta_outputs(&options, &annotation, genome.as_ref())?;
        sort_annotation(&mut annotation, &options.ref_sort_order, options.keep_genes)?;
    }
    let locus_rank = annotation
        .ref_order
        .iter()
        .enumerate()
        .map(|(index, seqid)| (seqid.as_str(), index))
        .collect::<std::collections::BTreeMap<_, _>>();
    loci.sort_by(|left, right| {
        (
            locus_rank
                .get(left.seqid.as_str())
                .copied()
                .unwrap_or(usize::MAX),
            left.start,
            left.end,
            &left.id,
        )
            .cmp(&(
                locus_rank
                    .get(right.seqid.as_str())
                    .copied()
                    .unwrap_or(usize::MAX),
                right.start,
                right.end,
                &right.id,
            ))
    });
    if options.expose_warnings {
        eprint!(
            "Command line was:\n{command_line}\n   .. loaded {} genomic features from {}\n",
            annotation.transcripts.len(),
            options.input.display()
        );
    }

    match options.main_output {
        MainOutput::Gff3 => {
            let gff3_options = gff3::Gff3Options {
                version: VERSION,
                command_line: &command_line,
                track_label: options.track_label.as_deref(),
                attrs_filter: options.attrs.as_deref(),
                keep_all_attrs: options.keep_all_attrs,
                gather_exon_attrs: options.gather_exon_attrs,
                keep_exon_attrs: options.keep_exon_attrs,
                keep_genes: options.keep_genes,
                keep_comments: options.keep_comments,
                decode_attrs: options.decode_attrs,
            };

            if let Some(output) = &options.output {
                let mut file = File::create(output).map_err(|_| {
                    CompatError::new(format!("Error creating file: {}\n", output.display()), 1)
                })?;

                gff3::write_gff3(&mut file, &annotation, &loci, &gff3_options)
                    .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
            } else if options.fasta_outputs.transcript.is_none()
                && options.fasta_outputs.unspliced.is_none()
                && options.fasta_outputs.cds.is_none()
                && options.fasta_outputs.protein.is_none()
            {
                let stdout = io::stdout();
                let mut handle = stdout.lock();
                gff3::write_gff3(&mut handle, &annotation, &loci, &gff3_options)
                    .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
            }
        }
        MainOutput::Gtf => {
            if let Some(output) = &options.output {
                let mut file = File::create(output).map_err(|_| {
                    CompatError::new(format!("Error creating file: {}\n", output.display()), 1)
                })?;

                gtf::write_gtf(&mut file, &annotation, options.track_label.as_deref())
                    .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
            } else {
                let stdout = io::stdout();
                let mut handle = stdout.lock();
                gtf::write_gtf(&mut handle, &annotation, options.track_label.as_deref())
                    .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
            }
        }
        MainOutput::Bed => {
            if let Some(output) = &options.output {
                let mut file = File::create(output).map_err(|_| {
                    CompatError::new(format!("Error creating file: {}\n", output.display()), 1)
                })?;
                bed::write_bed(&mut file, &annotation)
                    .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
            } else {
                let stdout = io::stdout();
                let mut handle = stdout.lock();
                bed::write_bed(&mut handle, &annotation)
                    .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
            }
        }
        MainOutput::Tlf => {
            if let Some(output) = &options.output {
                let mut file = File::create(output).map_err(|_| {
                    CompatError::new(format!("Error creating file: {}\n", output.display()), 1)
                })?;
                tlf::write_tlf(&mut file, &annotation, options.track_label.as_deref())
                    .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
            } else {
                let stdout = io::stdout();
                let mut handle = stdout.lock();
                tlf::write_tlf(&mut handle, &annotation, options.track_label.as_deref())
                    .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
            }
        }
        MainOutput::Table => {
            let format = options
                .table_format
                .as_deref()
                .ok_or_else(|| CompatError::new("Error: --table requires a format\n", 1))?;
            if let Some(output) = &options.output {
                let mut file = File::create(output).map_err(|_| {
                    CompatError::new(format!("Error creating file: {}\n", output.display()), 1)
                })?;

                table::write_table(&mut file, &annotation, format)
                    .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
            } else {
                let stdout = io::stdout();
                let mut handle = stdout.lock();
                table::write_table(&mut handle, &annotation, format)
                    .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
            }
        }
    }

    Ok(())
}

fn command_line(program: &str, args: &[String]) -> String {
    let display_program = if program.ends_with("gffread-rs") {
        workspace_root().join("gffread").display().to_string()
    } else {
        program.to_owned()
    };
    if args.is_empty() {
        display_program
    } else {
        format!("{display_program} {}", args.join(" "))
    }
}

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .canonicalize()
        .expect("workspace root must resolve")
}

fn write_fasta_outputs(
    options: &RuntimeOptions,
    annotation: &Annotation,
    genome: Option<&gffread_core::fasta::Genome>,
) -> Result<(), CompatError> {
    let Some(genome) = genome else {
        return Ok(());
    };

    if let Some(path) = &options.fasta_outputs.transcript {
        let file = File::create(path).map_err(|_| {
            CompatError::new(format!("Error creating file: {}\n", path.display()), 1)
        })?;
        let mut file = BufWriter::new(file);
        write_transcript_fasta(
            &mut file,
            annotation,
            genome,
            options.fasta_outputs.padding,
            options.fasta_outputs.suppress_transcript_cds,
            options.fasta_outputs.write_exon_segments,
        )?;
    }

    if let Some(path) = &options.fasta_outputs.unspliced {
        let file = File::create(path).map_err(|_| {
            CompatError::new(format!("Error creating file: {}\n", path.display()), 1)
        })?;
        let mut file = BufWriter::new(file);
        write_unspliced_fasta(&mut file, annotation, genome, options.fasta_outputs.padding)?;
    }

    if let Some(path) = &options.fasta_outputs.cds {
        let file = File::create(path).map_err(|_| {
            CompatError::new(format!("Error creating file: {}\n", path.display()), 1)
        })?;
        let mut file = BufWriter::new(file);
        write_cds_fasta(
            &mut file,
            annotation,
            genome,
            options.fasta_outputs.write_exon_segments,
        )?;
    }

    if let Some(path) = &options.fasta_outputs.protein {
        let file = File::create(path).map_err(|_| {
            CompatError::new(format!("Error creating file: {}\n", path.display()), 1)
        })?;
        let mut file = BufWriter::new(file);
        write_protein_fasta(
            &mut file,
            annotation,
            genome,
            options.fasta_outputs.write_protein_star_stop,
            options.fasta_outputs.write_exon_segments,
            options.fasta_outputs.cds.is_none(),
        )?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::Path;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::run_outputs;
    use gffread_core::options::{
        FastaOutputs, InputFormat, MainOutput, RefSortOrder, RuntimeOptions,
    };

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
            program: "./gffread".to_owned(),
            expose_warnings: false,
            output: Some(output.to_path_buf()),
            track_label: None,
            main_output: MainOutput::Gff3,
            table_format: None,
            attrs: None,
            keep_all_attrs: false,
            gather_exon_attrs: false,
            keep_exon_attrs: false,
            keep_genes: false,
            keep_comments: false,
            decode_attrs: false,
            ref_sort_order: RefSortOrder::Input,
            genome,
            fasta_outputs,
            range_filter: None,
            id_filter: None,
            input_format: InputFormat::Auto,
            min_length: None,
            max_intron: None,
            coding_only: false,
            noncoding_only: false,
            multi_exon_only: false,
            no_pseudo: false,
            cluster: gffread_core::options::ClusterOptions::default(),
            input: input.to_path_buf(),
            inputs: vec![input.to_path_buf()],
            original_args: Vec::new(),
        }
    }

    #[test]
    fn protein_fasta_option_can_run_with_gff3_output() {
        let fixtures = copied_example_fixtures("protein_fasta");
        let output = temp_output_path("protein_fasta");
        let protein_fasta = fixtures.root.join("transcripts_prot.fa");
        let fai_path = fixtures.root.join("genome.fa.fai");

        let result = run_outputs(gff3_options(
            &fixtures.annotation,
            &output,
            FastaOutputs {
                protein: Some(protein_fasta.clone()),
                ..FastaOutputs::default()
            },
            Some(fixtures.genome.clone()),
        ));

        result.expect("protein FASTA should be supported");
        assert!(output.exists(), "GFF3 output should still be written");
        assert!(protein_fasta.exists(), "protein FASTA should be written");
        assert!(
            fai_path.exists(),
            "temporary genome FASTA index should be written"
        );

        let _ = fs::remove_file(output);
        let _ = fs::remove_file(protein_fasta);
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
        let fai_path = fixtures.root.join("genome.fa.fai");

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
        assert!(
            transcript_fasta.exists(),
            "transcript FASTA should be written"
        );
        assert!(cds_fasta.exists(), "CDS FASTA should be written");
        assert!(
            fai_path.exists(),
            "temporary genome FASTA index should be written"
        );

        let _ = fs::remove_file(output);
        let _ = fs::remove_file(transcript_fasta);
        let _ = fs::remove_file(cds_fasta);
    }

    #[test]
    fn invalid_genome_path_with_main_output_and_transcript_fasta_does_not_create_main_output() {
        let fixtures = copied_example_fixtures("missing_genome_with_main_output");
        let output = temp_output_path("missing_genome_with_main_output");
        let transcript_fasta = fixtures.root.join("transcripts.fa");
        let missing_genome = fixtures.root.join("missing.fa");

        let result = run_outputs(gff3_options(
            &fixtures.annotation,
            &output,
            FastaOutputs {
                transcript: Some(transcript_fasta.clone()),
                ..FastaOutputs::default()
            },
            Some(missing_genome),
        ));

        let err = result.expect_err("missing genome should fail");
        assert_eq!(err.exit_code, 1);
        assert!(
            !output.exists(),
            "main output file must not be created before FASTA prerequisites are validated"
        );
        assert!(
            !transcript_fasta.exists(),
            "transcript FASTA should not be created when genome loading fails"
        );
    }
}
