use gffread_test_harness::CompatCase;

fn candidate() -> &'static str {
    env!("CARGO_BIN_EXE_gffread-rs")
}

#[test]
fn unspliced_fasta_matches_oracle() {
    CompatCase::new("unspliced_fasta")
        .in_examples()
        .args([
            "-u",
            "transcript_span.fa",
            "-g",
            "fasta_option_genome.fa",
            "fasta_option_case.gff",
        ])
        .expected_files(["transcript_span.fa"])
        .assert_matches_oracle(candidate())
        .expect("-u output must match oracle");
}

#[test]
fn unspliced_fasta_with_padding_matches_oracle() {
    CompatCase::new("unspliced_fasta_padding")
        .in_examples()
        .args([
            "-u",
            "transcript_span.fa",
            "--w-add",
            "2",
            "-g",
            "fasta_option_genome.fa",
            "fasta_option_case.gff",
        ])
        .expected_files(["transcript_span.fa"])
        .assert_matches_oracle(candidate())
        .expect("-u --w-add output must match oracle");
}

#[test]
fn transcript_fasta_without_cds_annotation_matches_oracle() {
    CompatCase::new("transcript_fasta_no_cds")
        .in_examples()
        .args([
            "-w",
            "transcripts.fa",
            "--w-nocds",
            "-g",
            "fasta_option_genome.fa",
            "fasta_option_case.gff",
        ])
        .expected_files(["transcripts.fa"])
        .assert_matches_oracle(candidate())
        .expect("-w --w-nocds output must match oracle");
}

#[test]
fn transcript_fasta_with_padding_matches_oracle() {
    CompatCase::new("transcript_fasta_padding")
        .in_examples()
        .args([
            "-w",
            "transcripts.fa",
            "--w-add",
            "2",
            "-g",
            "fasta_option_genome.fa",
            "fasta_option_case.gff",
        ])
        .expected_files(["transcripts.fa"])
        .assert_matches_oracle(candidate())
        .expect("-w --w-add output must match oracle");
}

#[test]
fn protein_fasta_with_star_stop_matches_oracle() {
    CompatCase::new("protein_fasta_star_stop")
        .in_examples()
        .args([
            "-y",
            "proteins.fa",
            "-S",
            "-g",
            "fasta_option_genome.fa",
            "fasta_option_case.gff",
        ])
        .expected_files(["proteins.fa"])
        .assert_matches_oracle(candidate())
        .expect("-y -S output must match oracle");
}

#[test]
fn transcript_fasta_with_projected_segments_matches_oracle() {
    CompatCase::new("transcript_fasta_with_segments")
        .in_examples()
        .args([
            "-W",
            "-w",
            "transcripts.fa",
            "-g",
            "fasta_option_genome.fa",
            "fasta_option_case.gff",
        ])
        .expected_files(["transcripts.fa"])
        .assert_matches_oracle(candidate())
        .expect("-W -w output must match oracle");
}

#[test]
fn transcript_fasta_with_segments_and_padding_matches_oracle() {
    CompatCase::new("transcript_fasta_with_segments_padding")
        .in_examples()
        .args([
            "-W",
            "-w",
            "transcripts.fa",
            "--w-add",
            "2",
            "-g",
            "fasta_option_genome.fa",
            "fasta_option_case.gff",
        ])
        .expected_files(["transcripts.fa"])
        .assert_matches_oracle(candidate())
        .expect("-W -w --w-add output must match oracle");
}

#[test]
fn protein_fasta_with_projected_segments_matches_oracle() {
    CompatCase::new("protein_fasta_with_segments")
        .in_examples()
        .args([
            "-W",
            "-y",
            "proteins.fa",
            "-g",
            "fasta_option_genome.fa",
            "fasta_option_case.gff",
        ])
        .expected_files(["proteins.fa"])
        .assert_matches_oracle(candidate())
        .expect("-W -y output must match oracle");
}
