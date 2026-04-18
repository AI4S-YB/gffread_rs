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
fn plain_gff3_stdout_matches_oracle() {
    CompatCase::new("example_gff3_stdout")
        .in_examples()
        .args(["annotation.gff"])
        .assert_matches_oracle(candidate())
        .expect("plain GFF3 stdout must match oracle");
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
fn missing_cds_phase_gff3_matches_oracle() {
    CompatCase::new("missing_cds_phase_gff3")
        .in_examples()
        .args(["phase_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("missing CDS phase must be recomputed like oracle");
}

#[test]
fn missing_cds_phase_gtf_matches_oracle() {
    CompatCase::new("missing_cds_phase_gtf")
        .in_examples()
        .args(["-T", "phase_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("missing CDS phase and GTF attributes must match oracle");
}

#[test]
fn edge_case_gff3_matches_oracle() {
    CompatCase::new("edge_case_gff3")
        .in_examples()
        .args(["edge_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("parent gene_name inheritance and exon merging must match oracle");
}

#[test]
fn edge_case_gtf_matches_oracle() {
    CompatCase::new("edge_case_gtf")
        .in_examples()
        .args(["-T", "edge_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("GTF exon phase and parent gene_name inheritance must match oracle");
}

#[test]
fn cds_extending_exon_gff3_matches_oracle() {
    CompatCase::new("cds_extending_exon_gff3")
        .in_examples()
        .args(["cds_extends_exon_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("exons should expand to include CDS bounds like oracle");
}

#[test]
fn cds_extending_exon_gtf_matches_oracle() {
    CompatCase::new("cds_extending_exon_gtf")
        .in_examples()
        .args(["-T", "cds_extends_exon_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("GTF exons should expand to include CDS bounds like oracle");
}

#[test]
fn exon_score_gff3_matches_oracle() {
    CompatCase::new("exon_score_gff3")
        .in_examples()
        .args(["segment_score_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("exon score column should match oracle");
}

#[test]
fn exon_score_gtf_matches_oracle() {
    CompatCase::new("exon_score_gtf")
        .in_examples()
        .args(["-T", "segment_score_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("GTF exon score column should match oracle");
}

#[test]
fn non_mrna_transcript_like_features_match_oracle() {
    CompatCase::new("non_mrna_transcript_like_features")
        .in_examples()
        .args(["non_mrna_transcript_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("non-mRNA transcript-like GFF3 features must match oracle");
}

#[test]
fn non_mrna_transcript_like_features_gtf_match_oracle() {
    CompatCase::new("non_mrna_transcript_like_features_gtf")
        .in_examples()
        .args(["-T", "non_mrna_transcript_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("non-mRNA transcript-like GTF conversion must match oracle");
}

#[test]
fn gene_name_equal_to_gene_id_gff3_matches_oracle() {
    CompatCase::new("gene_name_equal_to_gene_id_gff3")
        .in_examples()
        .args(["gene_name_id_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("gene_name matching gene_id should be suppressed like oracle");
}

#[test]
fn gene_name_equal_to_gene_id_gtf_matches_oracle() {
    CompatCase::new("gene_name_equal_to_gene_id_gtf")
        .in_examples()
        .args(["-T", "gene_name_id_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("GTF gene_name matching gene_id should be suppressed like oracle");
}

#[test]
fn duplicate_transcript_id_gff3_matches_oracle() {
    CompatCase::new("duplicate_transcript_id_gff3")
        .in_examples()
        .args(["duplicate_id_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("duplicate transcript IDs should attach children like oracle");
}

#[test]
fn duplicate_transcript_id_gtf_matches_oracle() {
    CompatCase::new("duplicate_transcript_id_gtf")
        .in_examples()
        .args(["-T", "duplicate_id_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("duplicate transcript IDs should attach children like oracle in GTF");
}

#[test]
fn duplicate_transcript_id_distinct_loci_gff3_matches_oracle() {
    CompatCase::new("duplicate_transcript_id_distinct_loci_gff3")
        .in_examples()
        .args(["duplicate_id_locus_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("duplicate transcript IDs at distinct loci should attach children like oracle");
}

#[test]
fn duplicate_transcript_id_distinct_loci_gtf_matches_oracle() {
    CompatCase::new("duplicate_transcript_id_distinct_loci_gtf")
        .in_examples()
        .args(["-T", "duplicate_id_locus_case.gff"])
        .assert_matches_oracle(candidate())
        .expect(
            "duplicate transcript IDs at distinct loci should attach children like oracle in GTF",
        );
}

#[test]
fn duplicate_transcript_id_same_ref_gff3_matches_oracle() {
    CompatCase::new("duplicate_transcript_id_same_ref_gff3")
        .in_examples()
        .args(["duplicate_id_same_ref_case.gff"])
        .assert_matches_oracle(candidate())
        .expect(
            "duplicate transcript IDs on the same reference should attach children like oracle",
        );
}

#[test]
fn duplicate_transcript_id_same_ref_gtf_matches_oracle() {
    CompatCase::new("duplicate_transcript_id_same_ref_gtf")
        .in_examples()
        .args(["-T", "duplicate_id_same_ref_case.gff"])
        .assert_matches_oracle(candidate())
        .expect(
            "duplicate transcript IDs on the same reference should attach children like oracle in GTF",
        );
}

#[test]
fn equal_key_duplicate_transcript_order_matches_oracle() {
    CompatCase::new("equal_key_duplicate_transcript_order")
        .in_examples()
        .args(["equal_key_duplicate_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("same-key duplicate transcripts must be ordered like gclib");
}

#[test]
fn duplicate_transcript_fasta_order_matches_oracle() {
    CompatCase::new("duplicate_transcript_fasta_order")
        .in_examples()
        .args([
            "-w",
            "transcripts.fa",
            "-g",
            "duplicate_fasta_order_genome.fa",
            "duplicate_fasta_order_case.gff",
        ])
        .expected_files(["transcripts.fa"])
        .assert_matches_oracle(candidate())
        .expect("FASTA output must follow the oracle loader order for duplicate transcript IDs");
}

#[test]
fn bed_conversion_matches_oracle() {
    CompatCase::new("example_bed")
        .in_examples()
        .args(["--bed", "annotation.gff"])
        .assert_matches_oracle(candidate())
        .expect("BED output must match oracle");
}

#[test]
fn bed_conversion_to_file_matches_oracle() {
    CompatCase::new("example_bed_file")
        .in_examples()
        .args(["--bed", "-o", "annotation.bed", "annotation.gff"])
        .expected_files(["annotation.bed"])
        .assert_matches_oracle(candidate())
        .expect("BED file output must match oracle");
}

#[test]
fn tlf_conversion_matches_oracle() {
    CompatCase::new("example_tlf")
        .in_examples()
        .args(["--tlf", "annotation.gff"])
        .assert_matches_oracle(candidate())
        .expect("TLF output must match oracle");
}

#[test]
fn tlf_conversion_to_file_matches_oracle() {
    CompatCase::new("example_tlf_file")
        .in_examples()
        .args(["--tlf", "-o", "annotation.tlf", "annotation.gff"])
        .expected_files(["annotation.tlf"])
        .assert_matches_oracle(candidate())
        .expect("TLF file output must match oracle");
}

#[test]
fn default_reference_order_matches_oracle() {
    CompatCase::new("sort_default")
        .in_examples()
        .args(["sort_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("default reference ordering must match oracle");
}

#[test]
fn sort_alpha_output_matches_oracle() {
    CompatCase::new("sort_alpha_output")
        .in_examples()
        .args(["--sort-alpha", "sort_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("--sort-alpha ordering must match oracle");
}

#[test]
fn sort_by_output_matches_oracle() {
    CompatCase::new("sort_by_output")
        .in_examples()
        .args(["--sort-by", "sort_refs.lst", "sort_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("--sort-by ordering must match oracle");
}

#[test]
fn track_label_gff3_matches_oracle() {
    CompatCase::new("track_label_gff3")
        .in_examples()
        .args(["-t", "custom", "annotation.gff"])
        .assert_matches_oracle(candidate())
        .expect("-t GFF3 output must match oracle");
}

#[test]
fn track_label_tlf_matches_oracle() {
    CompatCase::new("track_label_tlf")
        .in_examples()
        .args(["--tlf", "-t", "custom", "annotation.gff"])
        .assert_matches_oracle(candidate())
        .expect("-t TLF output must match oracle");
}

#[test]
fn decoded_full_attributes_match_oracle() {
    CompatCase::new("decode_full_attrs")
        .in_examples()
        .args(["-F", "-D", "decode_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("-F -D output must match oracle");
}

#[test]
fn full_attributes_reduce_redundant_exon_and_cds_attrs_matches_oracle() {
    CompatCase::new("reduce_exon_attrs")
        .in_examples()
        .args(["-F", "exon_attrs_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("-F output must match oracle");
}

#[test]
fn g_option_moves_first_exon_and_cds_attrs_to_transcript_matches_oracle() {
    CompatCase::new("move_exon_attrs_to_transcript")
        .in_examples()
        .args(["-G", "exon_attrs_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("-G output must match oracle");
}

#[test]
fn keep_exon_attrs_with_full_attributes_matches_oracle() {
    CompatCase::new("keep_exon_attrs")
        .in_examples()
        .args(["-F", "--keep-exon-attrs", "exon_attrs_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("-F --keep-exon-attrs output must match oracle");
}

#[test]
fn keep_genes_matches_oracle() {
    CompatCase::new("keep_genes")
        .in_examples()
        .args(["--keep-genes", "keep_genes_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("--keep-genes output must match oracle");
}

#[test]
fn keep_genes_with_full_attributes_matches_oracle() {
    CompatCase::new("keep_genes_full_attrs")
        .in_examples()
        .args(["--keep-genes", "-F", "keep_genes_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("--keep-genes -F output must match oracle");
}

#[test]
fn keep_comments_matches_oracle() {
    CompatCase::new("keep_comments")
        .in_examples()
        .args(["--keep-comments", "comments_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("--keep-comments output must match oracle");
}

#[test]
fn bed_input_extension_matches_oracle() {
    CompatCase::new("bed_input_extension")
        .in_examples()
        .args(["input_case.bed"])
        .assert_matches_oracle(candidate())
        .expect(".bed input parsing must match oracle");
}

#[test]
fn bed_input_flag_matches_oracle() {
    CompatCase::new("bed_input_flag")
        .in_examples()
        .args(["--in-bed", "input_case.bed"])
        .assert_matches_oracle(candidate())
        .expect("--in-bed parsing must match oracle");
}

#[test]
fn tlf_input_extension_matches_oracle() {
    CompatCase::new("tlf_input_extension")
        .in_examples()
        .args(["input_case.tlf"])
        .assert_matches_oracle(candidate())
        .expect(".tlf input parsing must match oracle");
}

#[test]
fn tlf_input_flag_matches_oracle() {
    CompatCase::new("tlf_input_flag")
        .in_examples()
        .args(["--in-tlf", "input_case.tlf"])
        .assert_matches_oracle(candidate())
        .expect("--in-tlf parsing must match oracle");
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

#[test]
fn protein_fasta_matches_oracle() {
    CompatCase::new("example_protein_fasta")
        .in_examples()
        .args([
            "-y",
            "transcripts_prot.fa",
            "-g",
            "genome.fa",
            "annotation.gff",
        ])
        .expected_files(["transcripts_prot.fa"])
        .assert_matches_oracle(candidate())
        .expect("protein FASTA must match oracle");
}

#[test]
fn combined_transcript_and_protein_fasta_match_oracle() {
    CompatCase::new("example_combined_fasta")
        .in_examples()
        .args([
            "-w",
            "transcripts.fa",
            "-y",
            "transcripts_prot.fa",
            "-g",
            "genome.fa",
            "annotation.gff",
        ])
        .expected_files(["transcripts.fa", "transcripts_prot.fa"])
        .assert_matches_oracle(candidate())
        .expect("combined transcript and protein FASTA must match oracle");
}
