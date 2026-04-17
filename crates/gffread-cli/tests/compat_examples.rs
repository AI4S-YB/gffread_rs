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

#[test]
fn transcript_fasta_matches_oracle() {
    CompatCase::new("example_transcript_fasta")
        .in_examples()
        .args(["-w", "transcripts.fa", "-g", "genome.fa", "annotation.gff"])
        .expected_files(["transcripts.fa"])
        .assert_matches_oracle(candidate())
        .expect("transcript FASTA must match oracle");
}

#[test]
fn cds_fasta_with_projected_segments_matches_oracle() {
    CompatCase::new("example_cds_fasta")
        .in_examples()
        .args([
            "-W",
            "-x",
            "transcripts_CDS.fa",
            "-g",
            "genome.fa",
            "annotation.gff",
        ])
        .expected_files(["transcripts_CDS.fa"])
        .assert_matches_oracle(candidate())
        .expect("CDS FASTA must match oracle");
}
