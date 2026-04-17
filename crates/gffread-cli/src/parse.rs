use std::fs::File;
use std::path::PathBuf;

use gffread_core::compat::CompatError;
use gffread_core::options::{FastaOutputs, MainOutput, RuntimeOptions};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CommandMode {
    Version,
    Help,
    Run(RuntimeOptions),
}

pub fn parse_args(args: Vec<String>) -> Result<CommandMode, CompatError> {
    if args == ["--version"] {
        return Ok(CommandMode::Version);
    }

    if args.is_empty() || args.iter().any(|arg| arg == "-h" || arg == "--help") {
        return Ok(CommandMode::Help);
    }

    let mut expose_warnings = false;
    let mut output = None;
    let mut main_output = MainOutput::Gff3;
    let mut table_format = None;
    let mut genome = None;
    let mut fasta_outputs = FastaOutputs::default();
    let mut sort_alpha = false;
    let mut sort_by = None::<String>;
    let mut inputs = Vec::<PathBuf>::new();

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-E" | "-v" => expose_warnings = true,
            "-T" | "--gtf" => main_output = MainOutput::Gtf,
            "-W" => fasta_outputs.write_exon_segments = true,
            "-o" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option -o requires an argument\n", 1)
                })?;
                output = Some(PathBuf::from(value));
            }
            "-g" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option -g requires an argument\n", 1)
                })?;
                genome = Some(PathBuf::from(value));
            }
            "-w" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option -w requires an argument\n", 1)
                })?;
                fasta_outputs.transcript = Some(PathBuf::from(value));
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
            "--table" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| {
                    CompatError::new("Error: option --table requires an argument\n", 1)
                })?;
                table_format = Some(value.clone());
                main_output = MainOutput::Table;
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

    if genome.is_none()
        && (fasta_outputs.transcript.is_some()
            || fasta_outputs.cds.is_some()
            || fasta_outputs.protein.is_some())
    {
        return Err(CompatError::new(
            "Error: -g option is required for options -w/x/y/u/V/N/M !\n",
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
        expose_warnings,
        output,
        main_output,
        table_format,
        genome,
        fasta_outputs,
        input,
        inputs,
        original_args: args,
    }))
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

        let mode = parse_args(args.clone()).expect("ordered positional inputs should parse");

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
            }
            other => panic!("expected run mode, got {other:?}"),
        }
    }
}
