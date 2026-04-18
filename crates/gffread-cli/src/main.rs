mod help;
mod parse;
mod process;

fn main() {
    let mut args = std::env::args();
    let program = args.next().unwrap_or_else(|| "gffread-rs".to_owned());
    let rest = args.collect::<Vec<_>>();
    std::process::exit(process::run_process(program, rest));
}
