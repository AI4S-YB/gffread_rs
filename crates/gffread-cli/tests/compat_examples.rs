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

#[test]
fn gtf_conversion_matches_oracle() {
    CompatCase::new("example_gtf")
        .in_examples()
        .args(["-T", "-o", "annotation.gtf", "annotation.gff"])
        .expected_files(["annotation.gtf"])
        .assert_matches_oracle(candidate())
        .expect("GTF output must match oracle");
}

#[test]
fn table_conversion_matches_oracle() {
    CompatCase::new("example_table")
        .in_examples()
        .args([
            "--table",
            "@id,@chr,@start,@end,@strand,@exons,Name,gene,product",
            "-o",
            "annotation.tbl",
            "annotation.gff",
        ])
        .expected_files(["annotation.tbl"])
        .assert_matches_oracle(candidate())
        .expect("table output must match oracle");
}
