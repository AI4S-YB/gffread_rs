use std::io::{self, Write};

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
        Ok(CommandMode::Run(_)) => {
            let _ =
                io::stderr().write_all(b"Error: Rust phase-one output path is not implemented\n");
            1
        }
        Err(err) => {
            let _ = io::stderr().write_all(err.message.as_bytes());
            err.exit_code
        }
    }
}
