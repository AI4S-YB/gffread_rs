use gffread_test_harness::CompatCase;

fn candidate() -> &'static str {
    env!("CARGO_BIN_EXE_gffread-rs")
}

#[test]
fn coding_only_filter_matches_oracle() {
    CompatCase::new("coding_only")
        .in_examples()
        .args(["-C", "coding_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("-C output must match oracle");
}

#[test]
fn noncoding_only_filter_matches_oracle() {
    CompatCase::new("noncoding_only")
        .in_examples()
        .args(["--nc", "coding_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("--nc output must match oracle");
}

#[test]
fn multi_exon_only_filter_matches_oracle() {
    CompatCase::new("multi_exon_only")
        .in_examples()
        .args(["-U", "coding_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("-U output must match oracle");
}

#[test]
fn advertised_add_hascds_error_matches_oracle() {
    CompatCase::new("add_hascds_error")
        .in_examples()
        .args(["--add-hasCDS", "coding_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("--add-hasCDS error surface must match oracle");
}
