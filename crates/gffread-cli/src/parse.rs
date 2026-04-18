use std::fs::File;
use std::path::PathBuf;

use gffread_core::compat::CompatError;
use gffread_core::options::{
    ClusterOptions, FastaOutputs, IdFilter, IdFilterMode, InputFormat, MainOutput, RangeFilter,
    RefSortOrder, RuntimeOptions,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CommandMode {
    Version,
    Help,
    UsageError(CompatError),
    Run(RuntimeOptions),
}

pub fn parse_args(program: String, args: Vec<String>) -> Result<CommandMode, CompatError> {
    if args == ["--version"] {
        return Ok(CommandMode::Version);
    }

    if args.is_empty() || args.iter().any(|arg| arg == "-h" || arg == "--help") {
        return Ok(CommandMode::Help);
    }

    let mut expose_warnings = false;
    let mut output = None;
    let mut track_label = None;
    let mut main_output = MainOutput::Gff3;
    let mut table_format = None;
    let mut attrs = None;
    let mut keep_all_attrs = false;
    let mut gather_exon_attrs = false;
    let mut keep_exon_attrs = false;
    let mut keep_genes = false;
    let mut keep_comments = false;
    let mut decode_attrs = false;
    let mut genome = None;
    let mut fasta_outputs = FastaOutputs::default();
    let mut sort_alpha = false;
    let mut sort_by = None::<String>;
    let mut range_filter = None;
    let mut range_within = false;
    let mut id_filter = None;
    let mut input_format = InputFormat::Auto;
    let mut min_length = None;
    let mut max_intron = None;
    let mut coding_only = false;
    let mut noncoding_only = false;
    let mut multi_exon_only = false;
    let mut no_pseudo = false;
    let mut cluster = ClusterOptions::default();
    let mut inputs = Vec::<PathBuf>::new();

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-E" | "-v" => expose_warnings = true,
            "-T" | "--gtf" => main_output = MainOutput::Gtf,
            "--bed" => main_output = MainOutput::Bed,
            "--tlf" => main_output = MainOutput::Tlf,
            "-D" => decode_attrs = true,
            "-F" => keep_all_attrs = true,
            "-G" => {
                gather_exon_attrs = true;
                keep_all_attrs = true;
            }
            "--keep-exon-attrs" => keep_exon_attrs = true,
            "--keep-genes" => keep_genes = true,
            "--keep-comments" => keep_comments = true,
            "--in-bed" => input_format = InputFormat::Bed,
            "--in-tlf" => input_format = InputFormat::Tlf,
            "-W" => fasta_outputs.write_exon_segments = true,
            "-R" => range_within = true,
            "-M" | "--merge" => cluster.merge = true,
            "-Q" => {
                cluster.merge = true;
                cluster.relax_boundary_containment = true;
            }
            "-K" => {
                cluster.merge = true;
                cluster.collapse_contained = true;
            }
            "--cluster-only" => {
                cluster.merge = true;
                cluster.cluster_only = true;
            }
            "-o" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option -o requires an argument\n", 1)
                })?;
                output = Some(PathBuf::from(value));
            }
            "-t" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option -t requires an argument\n", 1)
                })?;
                track_label = Some(value.clone());
            }
            "-r" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option -r requires an argument\n", 1)
                })?;
                range_filter = Some(parse_range(value)?);
            }
            "-C" => coding_only = true,
            "--nc" => noncoding_only = true,
            "-U" => multi_exon_only = true,
            "-g" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option -g requires an argument\n", 1)
                })?;
                genome = Some(PathBuf::from(value));
            }
            "-i" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option -i requires an argument\n", 1)
                })?;
                max_intron = Some(parse_u64_arg(value, "-i")?);
            }
            "-l" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option -l requires an argument\n", 1)
                })?;
                min_length = Some(parse_u64_arg(value, "-l")?);
            }
            "-w" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option -w requires an argument\n", 1)
                })?;
                fasta_outputs.transcript = Some(PathBuf::from(value));
            }
            "-u" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option -u requires an argument\n", 1)
                })?;
                fasta_outputs.unspliced = Some(PathBuf::from(value));
            }
            "-x" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option -x requires an argument\n", 1)
                })?;
                fasta_outputs.cds = Some(PathBuf::from(value));
            }
            "-y" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option -y requires an argument\n", 1)
                })?;
                fasta_outputs.protein = Some(PathBuf::from(value));
            }
            "-S" => fasta_outputs.write_protein_star_stop = true,
            "--table" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option --table requires an argument\n", 1)
                })?;
                table_format = Some(value.clone());
                main_output = MainOutput::Table;
            }
            "--attrs" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option --attrs requires an argument\n", 1)
                })?;
                attrs = Some(
                    value
                        .split(',')
                        .map(str::trim)
                        .filter(|part| !part.is_empty())
                        .map(str::to_owned)
                        .collect(),
                );
            }
            "--ids" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option --ids requires an argument\n", 1)
                })?;
                id_filter = Some(IdFilter {
                    path: PathBuf::from(value),
                    mode: IdFilterMode::Include,
                });
            }
            "--nids" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option --nids requires an argument\n", 1)
                })?;
                id_filter = Some(IdFilter {
                    path: PathBuf::from(value),
                    mode: IdFilterMode::Exclude,
                });
            }
            "--no-pseudo" => {
                no_pseudo = true;
                keep_all_attrs = true;
            }
            "--w-nocds" => fasta_outputs.suppress_transcript_cds = true,
            "--w-add" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option --w-add requires an argument\n", 1)
                })?;
                fasta_outputs.padding = parse_u64_arg(value, "--w-add")?;
            }
            "--add-hasCDS" => {
                return Ok(CommandMode::UsageError(CompatError::with_usage(
                    "Error: invalid argument '--add-hasCDS'\n",
                    1,
                )));
            }
            "--sort-alpha" => sort_alpha = true,
            "--sort-by" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option --sort-by requires an argument\n", 1)
                })?;
                sort_by = Some(value.clone());
            }
            value if value.starts_with('-') => {
                return Err(CompatError::new(
                    format!("Error: unknown option {value}\n"),
                    1,
                ));
            }
            value => inputs.push(PathBuf::from(value)),
        }

        i += 1;
    }

    if sort_alpha && sort_by.is_some() {
        return Err(CompatError::new(
            "Error: options --sort-by and --sort-alpha are mutually exclusive!\n",
            1,
        ));
    }

    if keep_exon_attrs && !keep_all_attrs {
        return Ok(CommandMode::UsageError(CompatError::new(
            "Error: option --keep-exon-attrs requires option -F !\n",
            0,
        )));
    }

    if range_within && range_filter.is_none() {
        return Err(CompatError::new("Error: option -R requires -r!\n", 1));
    }

    if let Some(range) = &mut range_filter {
        range.fully_within = range_within;
    }

    if genome.is_none()
        && (fasta_outputs.transcript.is_some()
            || fasta_outputs.unspliced.is_some()
            || fasta_outputs.cds.is_some()
            || fasta_outputs.protein.is_some())
    {
        return Err(CompatError::new(
            "Error: -g option is required for options -w/x/y/u/V/N/M !\n",
            1,
        ));
    }

    if fasta_outputs.padding > 0
        && fasta_outputs.transcript.is_none()
        && fasta_outputs.unspliced.is_none()
    {
        return Err(CompatError::new(
            "Error: --w-add option requires -w or -u option!\n",
            1,
        ));
    }

    if inputs.is_empty() {
        return Ok(CommandMode::Help);
    }

    for input in &inputs {
        File::open(input).map_err(|_| {
            CompatError::new(
                format!("Error: cannot open input file {}!\n", input.display()),
                1,
            )
        })?;
    }

    let input = inputs
        .first()
        .cloned()
        .expect("inputs must be non-empty after validation");

    Ok(CommandMode::Run(RuntimeOptions {
        program,
        expose_warnings,
        output,
        track_label,
        main_output,
        table_format,
        attrs,
        keep_all_attrs,
        gather_exon_attrs,
        keep_exon_attrs,
        keep_genes,
        keep_comments,
        decode_attrs,
        ref_sort_order: if let Some(path) = sort_by {
            RefSortOrder::List(PathBuf::from(path))
        } else if sort_alpha {
            RefSortOrder::Alpha
        } else {
            RefSortOrder::Input
        },
        genome,
        fasta_outputs,
        range_filter,
        id_filter,
        input_format,
        min_length,
        max_intron,
        coding_only,
        noncoding_only,
        multi_exon_only,
        no_pseudo,
        cluster,
        input,
        inputs,
        original_args: args,
    }))
}

fn parse_u64_arg(value: &str, option: &str) -> Result<u64, CompatError> {
    value
        .parse()
        .map_err(|_| CompatError::new(format!("Error: invalid value for option {option}\n"), 1))
}

fn parse_range(raw: &str) -> Result<RangeFilter, CompatError> {
    let (prefix, coords) = raw
        .split_once(':')
        .ok_or_else(|| CompatError::new("Error: invalid -r range format\n", 1))?;
    let (start, end) = coords
        .split_once('-')
        .ok_or_else(|| CompatError::new("Error: invalid -r range format\n", 1))?;

    let mut chars = prefix.chars();
    let first = chars.next().unwrap_or_default();
    let (strand, seqid) = if matches!(first, '+' | '-') {
        (Some(first), chars.as_str().to_owned())
    } else {
        (None, prefix.to_owned())
    };

    Ok(RangeFilter {
        seqid,
        strand,
        start: parse_u64_arg(start, "-r")?,
        end: parse_u64_arg(end, "-r")?,
        fully_within: false,
    })
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::{parse_args, CommandMode};
    use gffread_core::options::MainOutput;

    fn example_path(name: &str) -> String {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("..")
            .join("..")
            .join("examples")
            .join(name)
            .display()
            .to_string()
    }

    #[test]
    fn preserves_first_positional_input_as_primary() {
        let first_input = example_path("annotation.gff");
        let second_input = example_path("transcripts.gtf");
        let args = vec![first_input.clone(), second_input.clone()];

        let mode = parse_args("./gffread".to_owned(), args.clone())
            .expect("ordered positional inputs should parse");

        match mode {
            CommandMode::Run(options) => {
                assert_eq!(options.input, PathBuf::from(first_input.clone()));
                assert_eq!(
                    options.inputs,
                    vec![PathBuf::from(first_input), PathBuf::from(second_input)]
                );
                assert_eq!(options.main_output, MainOutput::Gff3);
                assert_eq!(options.output, None);
                assert_eq!(options.original_args, args);
                assert_eq!(options.program, "./gffread");
            }
            other => panic!("expected run mode, got {other:?}"),
        }
    }
}
