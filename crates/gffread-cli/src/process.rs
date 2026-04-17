use std::fs::File;
use std::io::{self, Write};
use std::path::PathBuf;

use gffread_core::compat::CompatError;
use gffread_core::emit::gff3;
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
        }
        Err(err) => {
            let _ = io::stderr().write_all(err.message.as_bytes());
            err.exit_code
        }
    }
}

fn run_outputs(options: RuntimeOptions) -> Result<(), CompatError> {
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
        MainOutput::Gtf | MainOutput::Table => Err(CompatError::new(
            "Error: selected output format is not implemented in this phase-one step\n",
            1,
        )),
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
