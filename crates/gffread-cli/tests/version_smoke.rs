use assert_cmd::Command;

#[test]
fn version_flag_prints_upstream_version() {
    Command::cargo_bin("gffread-rs")
        .unwrap()
        .arg("--version")
        .assert()
        .success()
        .stdout("0.12.9\n")
        .stderr("");
}
