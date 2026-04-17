use gffread_test_harness::CompatCase;

fn candidate() -> &'static str {
    env!("CARGO_BIN_EXE_gffread-rs")
}

#[test]
fn expose_to_simplified_gff3_matches_oracle() {
    CompatCase::new("example_gff3")
        .in_examples()
        .args(["-E", "-o", "ann_simple.gff", "annotation.gff"])
        .expected_files(["ann_simple.gff"])
        .assert_matches_oracle(candidate())
        .expect("simplified GFF3 output must match oracle");
}
