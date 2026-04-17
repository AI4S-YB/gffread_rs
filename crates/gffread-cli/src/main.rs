mod help;
mod parse;
mod process;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    std::process::exit(process::run_process(args));
}
