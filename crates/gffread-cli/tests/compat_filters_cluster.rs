use gffread_test_harness::CompatCase;

fn candidate() -> &'static str {
    env!("CARGO_BIN_EXE_gffread-rs")
}

#[test]
fn range_filter_matches_oracle() {
    CompatCase::new("range_filter")
        .in_examples()
        .args(["-r", "chr1:90-310", "filter_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("range filter output must match oracle");
}

#[test]
fn range_within_filter_matches_oracle() {
    CompatCase::new("range_within_filter")
        .in_examples()
        .args(["-r", "chr1:90-310", "-R", "filter_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("range-within filter output must match oracle");
}

#[test]
fn ids_filter_matches_oracle() {
    CompatCase::new("ids_filter")
        .in_examples()
        .args(["--ids", "filter_include.ids", "filter_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("include IDs filter output must match oracle");
}

#[test]
fn nids_filter_matches_oracle() {
    CompatCase::new("nids_filter")
        .in_examples()
        .args(["--nids", "filter_exclude.ids", "filter_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("exclude IDs filter output must match oracle");
}

#[test]
fn min_length_filter_matches_oracle() {
    CompatCase::new("min_length_filter")
        .in_examples()
        .args(["-l", "50", "filter_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("minimum length filter output must match oracle");
}

#[test]
fn max_intron_filter_matches_oracle() {
    CompatCase::new("max_intron_filter")
        .in_examples()
        .args(["-i", "30", "filter_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("maximum intron filter output must match oracle");
}

#[test]
fn no_pseudo_filter_matches_oracle() {
    CompatCase::new("no_pseudo_filter")
        .in_examples()
        .args(["--no-pseudo", "filter_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("no-pseudo filter output must match oracle");
}

#[test]
fn attrs_filter_matches_oracle() {
    CompatCase::new("attrs_filter")
        .in_examples()
        .args(["--attrs", "Name,product", "filter_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("attribute filter output must match oracle");
}

#[test]
fn merge_cluster_matches_oracle() {
    CompatCase::new("merge_cluster")
        .in_examples()
        .args(["-M", "cluster_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("merge clustering output must match oracle");
}

#[test]
fn merge_cluster_only_matches_oracle() {
    CompatCase::new("merge_cluster_only")
        .in_examples()
        .args(["-M", "--cluster-only", "cluster_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("cluster-only output must match oracle");
}

#[test]
fn merge_q_relaxes_boundary_containment_matches_oracle() {
    CompatCase::new("merge_q")
        .in_examples()
        .args(["-M", "-Q", "cluster_q_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("-M -Q output must match oracle");
}

#[test]
fn merge_k_collapses_contained_intron_chain_matches_oracle() {
    CompatCase::new("merge_k")
        .in_examples()
        .args(["-M", "-K", "cluster_k_case.gff"])
        .assert_matches_oracle(candidate())
        .expect("-M -K output must match oracle");
}
