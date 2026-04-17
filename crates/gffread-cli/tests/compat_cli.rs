use gffread_test_harness::CompatCase;

fn candidate() -> &'static str {
    env!("CARGO_BIN_EXE_gffread-rs")
}

#[test]
fn short_help_matches_oracle() {
    CompatCase::new("short_help")
        .args(["-h"])
        .assert_matches_oracle(candidate())
        .expect("-h must match oracle");
}

#[test]
fn long_help_matches_oracle() {
    CompatCase::new("long_help")
        .args(["--help"])
        .assert_matches_oracle(candidate())
        .expect("--help must match oracle");
}

#[test]
fn missing_genome_for_w_matches_oracle() {
    CompatCase::new("missing_genome_for_w")
        .in_examples()
        .args(["-w", "out.fa", "annotation.gff"])
        .assert_matches_oracle(candidate())
        .expect("missing -g error must match oracle");
}

#[test]
fn mutually_exclusive_sort_options_match_oracle() {
    CompatCase::new("sort_alpha_sort_by")
        .in_examples()
        .args(["--sort-alpha", "--sort-by", "refs.lst", "annotation.gff"])
        .assert_matches_oracle(candidate())
        .expect("sort option error must match oracle");
}

#[test]
fn missing_input_file_matches_oracle() {
    CompatCase::new("missing_input_file")
        .in_examples()
        .args(["does-not-exist.gff"])
        .assert_matches_oracle(candidate())
        .expect("missing input error must match oracle");
}
