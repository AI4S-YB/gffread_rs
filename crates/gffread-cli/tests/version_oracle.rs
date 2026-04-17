use gffread_test_harness::CompatCase;

#[test]
fn version_matches_cpp_oracle() {
    CompatCase::new("version")
        .args(["--version"])
        .assert_matches_oracle(env!("CARGO_BIN_EXE_gffread-rs"))
        .expect("version output must match oracle");
}
