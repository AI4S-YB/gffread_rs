fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();

    if args == ["--version"] {
        println!("{}", gffread_core::VERSION);
        return;
    }

    std::process::exit(1);
}
