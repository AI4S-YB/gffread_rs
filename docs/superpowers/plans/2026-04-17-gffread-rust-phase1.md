# GffRead Rust Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Rust `gffread-rs` executable and oracle-backed test harness that match the current local C++ `gffread v0.12.9` for phase-one CLI and example-output cases byte-for-byte.

**Architecture:** Add a Cargo workspace inside the existing `gffread/` repository. `gffread-rs` owns process compatibility and delegates parsing/output to `gffread-core`; `gffread-test-harness` builds/runs the local C++ oracle and compares exit code, stdout, stderr, and generated files.

**Tech Stack:** Rust 2021 edition, Cargo workspace, manual argv parsing with the standard library, `assert_cmd`, `tempfile`, `similar`, `thiserror`, and the existing C++ `make` build.

---

## File Structure

- Create: `.gitignore`
- Create: `Cargo.toml`
- Create: `crates/gffread-cli/Cargo.toml`
- Create: `crates/gffread-cli/src/main.rs`
- Create: `crates/gffread-cli/src/help.rs`
- Create: `crates/gffread-cli/src/parse.rs`
- Create: `crates/gffread-cli/src/process.rs`
- Create: `crates/gffread-core/Cargo.toml`
- Create: `crates/gffread-core/src/lib.rs`
- Create: `crates/gffread-core/src/compat.rs`
- Create: `crates/gffread-core/src/options.rs`
- Create: `crates/gffread-core/src/model.rs`
- Create: `crates/gffread-core/src/loader/mod.rs`
- Create: `crates/gffread-core/src/loader/gff.rs`
- Create: `crates/gffread-core/src/emit/mod.rs`
- Create: `crates/gffread-core/src/emit/gff3.rs`
- Create: `crates/gffread-core/src/emit/gtf.rs`
- Create: `crates/gffread-core/src/emit/table.rs`
- Create: `crates/gffread-core/src/fasta.rs`
- Create: `crates/gffread-test-harness/Cargo.toml`
- Create: `crates/gffread-test-harness/src/lib.rs`
- Create: `crates/gffread-cli/tests/version_smoke.rs`
- Create: `crates/gffread-cli/tests/version_oracle.rs`
- Create: `crates/gffread-cli/tests/compat_cli.rs`
- Create: `crates/gffread-cli/tests/compat_examples.rs`
- Create: `docs/compat/function-map.md`
- Modify: `Makefile`

The `gffread/examples` directory remains the canonical fixture source for phase-one end-to-end tests.

## Task 1: Bootstrap Workspace and Version Smoke Test

**Files:**

- Create: `.gitignore`
- Create: `Cargo.toml`
- Create: `crates/gffread-cli/Cargo.toml`
- Create: `crates/gffread-cli/src/main.rs`
- Create: `crates/gffread-core/Cargo.toml`
- Create: `crates/gffread-core/src/lib.rs`
- Create: `crates/gffread-test-harness/Cargo.toml`
- Create: `crates/gffread-cli/tests/version_smoke.rs`

- [ ] **Step 1: Write the failing workspace and version smoke test**

Create the initial workspace manifests and a binary that intentionally does not satisfy `--version`.

```gitignore
# .gitignore
/target/
```

```toml
# Cargo.toml
[workspace]
members = [
  "crates/gffread-cli",
  "crates/gffread-core",
  "crates/gffread-test-harness",
]
resolver = "2"

[workspace.package]
edition = "2021"
license = "MIT"
version = "0.1.0"

[workspace.dependencies]
assert_cmd = "2.0"
similar = "2.6"
tempfile = "3.15"
thiserror = "2.0"
```

```toml
# crates/gffread-core/Cargo.toml
[package]
name = "gffread-core"
edition.workspace = true
license.workspace = true
version.workspace = true

[dependencies]
thiserror.workspace = true
```

```rust
// crates/gffread-core/src/lib.rs
pub const VERSION: &str = "0.12.9";
```

```toml
# crates/gffread-cli/Cargo.toml
[package]
name = "gffread-rs"
edition.workspace = true
license.workspace = true
version.workspace = true

[dependencies]
gffread-core = { path = "../gffread-core" }

[dev-dependencies]
assert_cmd.workspace = true
```

```rust
// crates/gffread-cli/src/main.rs
fn main() {
    std::process::exit(1);
}
```

```toml
# crates/gffread-test-harness/Cargo.toml
[package]
name = "gffread-test-harness"
edition.workspace = true
license.workspace = true
version.workspace = true
```

```rust
// crates/gffread-cli/tests/version_smoke.rs
use assert_cmd::Command;

#[test]
fn version_flag_prints_upstream_version() {
    let mut cmd = Command::cargo_bin("gffread-rs").expect("gffread-rs binary");
    cmd.arg("--version")
        .assert()
        .success()
        .stdout("0.12.9\n")
        .stderr("");
}
```

- [ ] **Step 2: Run the smoke test and verify it fails**

Run: `cargo test -p gffread-rs version_flag_prints_upstream_version -- --exact`

Expected: FAIL because `gffread-rs` exits with status `1`.

- [ ] **Step 3: Implement the minimal `--version` path**

Replace `crates/gffread-cli/src/main.rs` with:

```rust
use gffread_core::VERSION;

fn main() {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args == ["--version"] {
        print!("{VERSION}\n");
        return;
    }

    std::process::exit(1);
}
```

- [ ] **Step 4: Run the smoke test and verify it passes**

Run: `cargo test -p gffread-rs version_flag_prints_upstream_version -- --exact`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add .gitignore Cargo.toml crates/gffread-cli crates/gffread-core crates/gffread-test-harness
git commit -m "feat: bootstrap rust workspace"
```

## Task 2: Build the Oracle Harness

**Files:**

- Create: `crates/gffread-test-harness/src/lib.rs`
- Create: `crates/gffread-cli/tests/version_oracle.rs`
- Modify: `crates/gffread-cli/Cargo.toml`
- Modify: `crates/gffread-test-harness/Cargo.toml`

- [ ] **Step 1: Write the failing oracle comparison test**

Update the harness manifest:

```toml
# crates/gffread-test-harness/Cargo.toml
[package]
name = "gffread-test-harness"
edition.workspace = true
license.workspace = true
version.workspace = true

[dependencies]
similar.workspace = true
tempfile.workspace = true
```

Update the CLI manifest:

```toml
# crates/gffread-cli/Cargo.toml
[package]
name = "gffread-rs"
edition.workspace = true
license.workspace = true
version.workspace = true

[dependencies]
gffread-core = { path = "../gffread-core" }

[dev-dependencies]
assert_cmd.workspace = true
gffread-test-harness = { path = "../gffread-test-harness" }
```

Create the test:

```rust
// crates/gffread-cli/tests/version_oracle.rs
use gffread_test_harness::CompatCase;

#[test]
fn version_matches_cpp_oracle() {
    CompatCase::new("version")
        .args(["--version"])
        .assert_matches_oracle(env!("CARGO_BIN_EXE_gffread-rs"))
        .expect("version output must match oracle");
}
```

- [ ] **Step 2: Run the oracle test and verify it fails to compile**

Run: `cargo test -p gffread-rs version_matches_cpp_oracle -- --exact`

Expected: FAIL with an unresolved import for `gffread_test_harness::CompatCase`.

- [ ] **Step 3: Implement `CompatCase` with process and file capture**

Create `crates/gffread-test-harness/src/lib.rs`:

```rust
use similar::{ChangeTag, TextDiff};
use std::collections::BTreeMap;
use std::error::Error;
use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use tempfile::TempDir;

#[derive(Debug, Clone)]
pub struct CompatCase {
    name: &'static str,
    args: Vec<String>,
    cwd: CaseCwd,
    expected_files: Vec<String>,
}

#[derive(Debug, Clone, Copy)]
enum CaseCwd {
    RepoRoot,
    Examples,
}

#[derive(Debug)]
struct CapturedRun {
    output: Output,
    files: BTreeMap<String, Vec<u8>>,
}

impl CompatCase {
    pub fn new(name: &'static str) -> Self {
        Self {
            name,
            args: Vec::new(),
            cwd: CaseCwd::RepoRoot,
            expected_files: Vec::new(),
        }
    }

    pub fn args<I, S>(mut self, args: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        self.args = args.into_iter().map(|arg| arg.as_ref().to_owned()).collect();
        self
    }

    pub fn in_examples(mut self) -> Self {
        self.cwd = CaseCwd::Examples;
        self
    }

    pub fn expected_files<I, S>(mut self, files: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        self.expected_files = files
            .into_iter()
            .map(|file| file.as_ref().to_owned())
            .collect();
        self
    }

    pub fn assert_matches_oracle<P>(self, candidate: P) -> Result<(), Box<dyn Error>>
    where
        P: AsRef<Path>,
    {
        let oracle = build_cpp_oracle()?;
        let oracle_run = self.run(&oracle)?;
        let candidate_run = self.run(candidate.as_ref())?;

        compare_status(self.name, &oracle_run.output, &candidate_run.output)?;
        compare_stream(self.name, "stdout", &oracle_run.output.stdout, &candidate_run.output.stdout)?;
        compare_stream(self.name, "stderr", &oracle_run.output.stderr, &candidate_run.output.stderr)?;
        compare_files(self.name, &oracle_run.files, &candidate_run.files)?;
        Ok(())
    }

    fn run(&self, binary: &Path) -> Result<CapturedRun, Box<dyn Error>> {
        let temp = TempDir::new()?;
        let cwd = match self.cwd {
            CaseCwd::RepoRoot => temp.path().to_path_buf(),
            CaseCwd::Examples => {
                copy_examples(temp.path())?;
                temp.path().join("examples")
            }
        };

        let output = Command::new(binary).args(&self.args).current_dir(&cwd).output()?;
        let mut files = BTreeMap::new();
        for expected in &self.expected_files {
            let path = cwd.join(expected);
            let bytes = fs::read(&path).map_err(|err| {
                format!(
                    "case {} expected output file {} was not readable: {}",
                    self.name, expected, err
                )
            })?;
            files.insert(expected.clone(), bytes);
        }
        Ok(CapturedRun { output, files })
    }
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(2)
        .expect("harness crate must live under crates/gffread-test-harness")
        .to_path_buf()
}

fn build_cpp_oracle() -> Result<PathBuf, Box<dyn Error>> {
    let root = repo_root();
    let oracle = root.join("gffread");
    if oracle.exists() {
        return Ok(oracle);
    }

    let status = Command::new("make").arg("-C").arg(&root).status()?;
    if !status.success() {
        return Err(format!("failed to build C++ oracle with make -C {}", root.display()).into());
    }
    Ok(oracle)
}

fn copy_examples(temp_root: &Path) -> Result<(), Box<dyn Error>> {
    let source = repo_root().join("examples");
    let dest = temp_root.join("examples");
    fs::create_dir_all(&dest)?;
    for entry in fs::read_dir(source)? {
        let entry = entry?;
        let file_type = entry.file_type()?;
        if file_type.is_file() {
            fs::copy(entry.path(), dest.join(entry.file_name()))?;
        }
    }
    Ok(())
}

fn compare_status(case: &str, oracle: &Output, candidate: &Output) -> Result<(), Box<dyn Error>> {
    if oracle.status.code() == candidate.status.code() {
        return Ok(());
    }
    Err(format!(
        "case {case} exit code differs: oracle {:?}, candidate {:?}",
        oracle.status.code(),
        candidate.status.code()
    )
    .into())
}

fn compare_stream(
    case: &str,
    stream: &str,
    oracle: &[u8],
    candidate: &[u8],
) -> Result<(), Box<dyn Error>> {
    if oracle == candidate {
        return Ok(());
    }
    Err(format!(
        "case {case} {stream} differs\n{}",
        unified_diff(stream, oracle, candidate)
    )
    .into())
}

fn compare_files(
    case: &str,
    oracle: &BTreeMap<String, Vec<u8>>,
    candidate: &BTreeMap<String, Vec<u8>>,
) -> Result<(), Box<dyn Error>> {
    if oracle.keys().ne(candidate.keys()) {
        return Err(format!(
            "case {case} output file set differs: oracle {:?}, candidate {:?}",
            oracle.keys().collect::<Vec<_>>(),
            candidate.keys().collect::<Vec<_>>()
        )
        .into());
    }
    for (name, oracle_bytes) in oracle {
        let candidate_bytes = candidate.get(name).expect("keys checked above");
        if oracle_bytes != candidate_bytes {
            return Err(format!(
                "case {case} file {name} differs\n{}",
                unified_diff(name, oracle_bytes, candidate_bytes)
            )
            .into());
        }
    }
    Ok(())
}

fn unified_diff(label: &str, oracle: &[u8], candidate: &[u8]) -> String {
    let oracle_text = String::from_utf8_lossy(oracle);
    let candidate_text = String::from_utf8_lossy(candidate);
    let diff = TextDiff::from_lines(&oracle_text, &candidate_text);
    let mut out = format!("--- oracle/{label}\n+++ candidate/{label}\n");
    for change in diff.iter_all_changes() {
        let sign = match change.tag() {
            ChangeTag::Delete => "-",
            ChangeTag::Insert => "+",
            ChangeTag::Equal => " ",
        };
        out.push_str(sign);
        out.push_str(change.value());
    }
    out
}
```

- [ ] **Step 4: Run the oracle test and verify it passes**

Run: `cargo test -p gffread-rs version_matches_cpp_oracle -- --exact`

Expected: PASS. If `../gclib` is absent, the first run may clone/build it through the existing `Makefile`.

- [ ] **Step 5: Commit**

```bash
git add crates/gffread-cli/Cargo.toml crates/gffread-cli/tests/version_oracle.rs crates/gffread-test-harness
git commit -m "test: add cpp oracle compatibility harness"
```

## Task 3: Reproduce CLI Help and Validation Surface

**Files:**

- Create: `crates/gffread-cli/src/help.rs`
- Create: `crates/gffread-cli/src/parse.rs`
- Create: `crates/gffread-cli/src/process.rs`
- Create: `crates/gffread-core/src/options.rs`
- Create: `crates/gffread-core/src/compat.rs`
- Modify: `crates/gffread-core/src/lib.rs`
- Modify: `crates/gffread-cli/src/main.rs`
- Create: `crates/gffread-cli/tests/compat_cli.rs`

- [ ] **Step 1: Write failing oracle tests for CLI surface cases**

Create `crates/gffread-cli/tests/compat_cli.rs`:

```rust
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
```

- [ ] **Step 2: Run the CLI tests and verify they fail**

Run: `cargo test -p gffread-rs --test compat_cli`

Expected: FAIL because `gffread-rs` still only handles `--version`.

- [ ] **Step 3: Add options and compatibility error types**

Create `crates/gffread-core/src/compat.rs`:

```rust
use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompatError {
    pub message: String,
    pub exit_code: i32,
}

impl CompatError {
    pub fn new(message: impl Into<String>, exit_code: i32) -> Self {
        Self {
            message: message.into(),
            exit_code,
        }
    }
}

impl fmt::Display for CompatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for CompatError {}
```

Create `crates/gffread-core/src/options.rs`:

```rust
use std::path::PathBuf;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MainOutput {
    Gff3,
    Gtf,
    Table,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FastaOutputs {
    pub transcript: Option<PathBuf>,
    pub cds: Option<PathBuf>,
    pub protein: Option<PathBuf>,
    pub write_exon_segments: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RuntimeOptions {
    pub expose_warnings: bool,
    pub output: Option<PathBuf>,
    pub main_output: MainOutput,
    pub table_format: Option<String>,
    pub genome: Option<PathBuf>,
    pub fasta_outputs: FastaOutputs,
    pub input: PathBuf,
    pub original_args: Vec<String>,
}
```

Update `crates/gffread-core/src/lib.rs`:

```rust
pub mod compat;
pub mod options;

pub const VERSION: &str = "0.12.9";
```

- [ ] **Step 4: Implement manual phase-one parsing and process dispatch**

Create `crates/gffread-cli/src/help.rs` with the exact usage bytes from the local C++ oracle:

```rust
pub const USAGE: &str = include_str!("usage.txt");
```

Create `crates/gffread-cli/src/parse.rs`:

```rust
use gffread_core::compat::CompatError;
use gffread_core::options::{FastaOutputs, MainOutput, RuntimeOptions};
use std::path::PathBuf;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CommandMode {
    Version,
    Help,
    Run(RuntimeOptions),
}

pub fn parse_args(args: Vec<String>) -> Result<CommandMode, CompatError> {
    if args == ["--version"] {
        return Ok(CommandMode::Version);
    }
    if args.iter().any(|arg| arg == "-h" || arg == "--help") || args.is_empty() {
        return Ok(CommandMode::Help);
    }

    let mut expose_warnings = false;
    let mut output = None;
    let mut main_output = MainOutput::Gff3;
    let mut table_format = None;
    let mut genome = None;
    let mut fasta_outputs = FastaOutputs::default();
    let mut sort_alpha = false;
    let mut sort_by = None::<String>;
    let mut input = None::<PathBuf>;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-E" | "-v" => expose_warnings = true,
            "-T" | "--gtf" => main_output = MainOutput::Gtf,
            "-W" => fasta_outputs.write_exon_segments = true,
            "-o" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| CompatError::new("Error: option -o requires an argument\n", 1))?;
                output = Some(PathBuf::from(value));
            }
            "-g" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| CompatError::new("Error: option -g requires an argument\n", 1))?;
                genome = Some(PathBuf::from(value));
            }
            "-w" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| CompatError::new("Error: option -w requires an argument\n", 1))?;
                fasta_outputs.transcript = Some(PathBuf::from(value));
            }
            "-x" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| CompatError::new("Error: option -x requires an argument\n", 1))?;
                fasta_outputs.cds = Some(PathBuf::from(value));
            }
            "-y" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| CompatError::new("Error: option -y requires an argument\n", 1))?;
                fasta_outputs.protein = Some(PathBuf::from(value));
            }
            "--table" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| CompatError::new("Error: option --table requires an argument\n", 1))?;
                table_format = Some(value.clone());
                main_output = MainOutput::Table;
            }
            "--sort-alpha" => sort_alpha = true,
            "--sort-by" => {
                i += 1;
                let value = args.get(i).ok_or_else(|| CompatError::new("Error: option --sort-by requires an argument\n", 1))?;
                sort_by = Some(value.clone());
            }
            value if value.starts_with('-') => {
                return Err(CompatError::new(format!("Error: unknown option {value}\n"), 1));
            }
            value => input = Some(PathBuf::from(value)),
        }
        i += 1;
    }

    if sort_alpha && sort_by.is_some() {
        return Err(CompatError::new(
            "Error: options --sort-by and --sort-alpha are mutually exclusive!\n",
            1,
        ));
    }

    if genome.is_none()
        && (fasta_outputs.transcript.is_some()
            || fasta_outputs.cds.is_some()
            || fasta_outputs.protein.is_some())
    {
        return Err(CompatError::new(
            "Error: -g option is required for options -w/x/y/u/V/N/M !\n",
            1,
        ));
    }

    let input = input.ok_or_else(|| CompatError::new(super::help::USAGE.to_owned(), 1))?;
    if !input.exists() {
        return Err(CompatError::new(
            format!("Error: cannot open input file {}!\n", input.display()),
            1,
        ));
    }

    Ok(CommandMode::Run(RuntimeOptions {
        expose_warnings,
        output,
        main_output,
        table_format,
        genome,
        fasta_outputs,
        input,
        original_args: args,
    }))
}
```

Create `crates/gffread-cli/src/process.rs`:

```rust
use crate::help::USAGE;
use crate::parse::{parse_args, CommandMode};
use gffread_core::VERSION;
use std::io::{self, Write};

pub fn run_process(args: Vec<String>) -> i32 {
    match parse_args(args) {
        Ok(CommandMode::Version) => {
            print!("{VERSION}\n");
            0
        }
        Ok(CommandMode::Help) => {
            let _ = io::stderr().write_all(USAGE.as_bytes());
            1
        }
        Ok(CommandMode::Run(_options)) => {
            eprint!("Error: Rust phase-one output path is not implemented\n");
            1
        }
        Err(err) => {
            let _ = io::stderr().write_all(err.message.as_bytes());
            err.exit_code
        }
    }
}
```

Replace `crates/gffread-cli/src/main.rs`:

```rust
mod help;
mod parse;
mod process;

fn main() {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    std::process::exit(process::run_process(args));
}
```

Create `crates/gffread-cli/src/usage.txt` by copying the exact C++ `USAGE` text from `gffread.cpp`. Confirm with the oracle test rather than visual inspection.

- [ ] **Step 5: Run CLI tests and adjust only byte differences**

Run: `cargo test -p gffread-rs --test compat_cli`

Expected: PASS after `usage.txt` and validation text match the C++ oracle. If a failure is only a wording or newline mismatch, change the compatibility renderer instead of changing the test.

- [ ] **Step 6: Commit**

```bash
git add crates/gffread-cli crates/gffread-core
git commit -m "feat: match phase one cli surface"
```

## Task 4: Parse Example GFF and Emit Simplified GFF3

**Files:**

- Create: `crates/gffread-core/src/model.rs`
- Create: `crates/gffread-core/src/loader/mod.rs`
- Create: `crates/gffread-core/src/loader/gff.rs`
- Create: `crates/gffread-core/src/emit/mod.rs`
- Create: `crates/gffread-core/src/emit/gff3.rs`
- Modify: `crates/gffread-core/src/lib.rs`
- Modify: `crates/gffread-cli/src/process.rs`
- Create: `crates/gffread-cli/tests/compat_examples.rs`

- [ ] **Step 1: Write the failing `-E` example parity test**

Create `crates/gffread-cli/tests/compat_examples.rs`:

```rust
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
```

- [ ] **Step 2: Run the `-E` example test and verify it fails**

Run: `cargo test -p gffread-rs expose_to_simplified_gff3_matches_oracle -- --exact`

Expected: FAIL because the run path still emits "not implemented".

- [ ] **Step 3: Add the phase-one model and GFF loader**

Create `crates/gffread-core/src/model.rs`:

```rust
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Annotation {
    pub transcripts: Vec<Transcript>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Transcript {
    pub seqid: String,
    pub source: String,
    pub feature: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub id: String,
    pub gene_id: Option<String>,
    pub gene_name: Option<String>,
    pub attrs: BTreeMap<String, String>,
    pub exons: Vec<Segment>,
    pub cds: Vec<Segment>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Segment {
    pub start: u64,
    pub end: u64,
    pub phase: String,
}

impl Transcript {
    pub fn exon_list(&self) -> String {
        self.exons
            .iter()
            .map(|seg| format!("{}-{}", seg.start, seg.end))
            .collect::<Vec<_>>()
            .join(",")
    }
}
```

Create `crates/gffread-core/src/loader/mod.rs`:

```rust
pub mod gff;
```

Create `crates/gffread-core/src/loader/gff.rs`:

```rust
use crate::compat::CompatError;
use crate::model::{Annotation, Segment, Transcript};
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

pub fn load_annotation(path: &Path) -> Result<Annotation, CompatError> {
    let text = fs::read_to_string(path)
        .map_err(|_| CompatError::new(format!("Error: cannot open input file {}!\n", path.display()), 1))?;
    parse_annotation(&text)
}

pub fn parse_annotation(text: &str) -> Result<Annotation, CompatError> {
    let mut transcripts: Vec<Transcript> = Vec::new();
    let mut by_id: BTreeMap<String, usize> = BTreeMap::new();

    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != 9 {
            continue;
        }
        let attrs = parse_attrs(fields[8]);
        let feature = fields[2];
        match feature {
            "mRNA" | "transcript" => {
                let id = attr(&attrs, "ID")
                    .or_else(|| attr(&attrs, "transcript_id"))
                    .ok_or_else(|| CompatError::new("Error: transcript record without ID\n", 1))?;
                let transcript = Transcript {
                    seqid: fields[0].to_owned(),
                    source: fields[1].to_owned(),
                    feature: feature.to_owned(),
                    start: fields[3].parse().map_err(|_| CompatError::new("Error: invalid transcript start\n", 1))?,
                    end: fields[4].parse().map_err(|_| CompatError::new("Error: invalid transcript end\n", 1))?,
                    strand: fields[6].chars().next().unwrap_or('.'),
                    id: id.clone(),
                    gene_id: attr(&attrs, "gene").or_else(|| attr(&attrs, "gene_id")),
                    gene_name: attr(&attrs, "gene_name"),
                    attrs,
                    exons: Vec::new(),
                    cds: Vec::new(),
                };
                by_id.insert(id, transcripts.len());
                transcripts.push(transcript);
            }
            "exon" | "CDS" => {
                let parent = attr(&attrs, "Parent")
                    .or_else(|| attr(&attrs, "transcript_id"))
                    .ok_or_else(|| CompatError::new("Error: child feature without Parent\n", 1))?;
                if let Some(index) = by_id.get(&parent).copied() {
                    let segment = Segment {
                        start: fields[3].parse().map_err(|_| CompatError::new("Error: invalid segment start\n", 1))?,
                        end: fields[4].parse().map_err(|_| CompatError::new("Error: invalid segment end\n", 1))?,
                        phase: fields[7].to_owned(),
                    };
                    if feature == "exon" {
                        transcripts[index].exons.push(segment);
                    } else {
                        transcripts[index].cds.push(segment);
                    }
                }
            }
            _ => {}
        }
    }

    Ok(Annotation { transcripts })
}

fn parse_attrs(raw: &str) -> BTreeMap<String, String> {
    let mut attrs = BTreeMap::new();
    for part in raw.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        if let Some((key, value)) = part.split_once('=') {
            attrs.insert(key.trim().to_owned(), value.trim().trim_matches('"').to_owned());
        } else if let Some((key, value)) = part.split_once(' ') {
            attrs.insert(key.trim().to_owned(), value.trim().trim_matches('"').to_owned());
        }
    }
    attrs
}

fn attr(attrs: &BTreeMap<String, String>, name: &str) -> Option<String> {
    attrs.get(name).filter(|value| !value.is_empty()).cloned()
}
```

Update `crates/gffread-core/src/lib.rs`:

```rust
pub mod compat;
pub mod loader;
pub mod model;
pub mod options;

pub const VERSION: &str = "0.12.9";
```

- [ ] **Step 4: Add GFF3 emitter and wire `-E -o` run path**

Create `crates/gffread-core/src/emit/mod.rs`:

```rust
pub mod gff3;
```

Create `crates/gffread-core/src/emit/gff3.rs`:

```rust
use crate::model::{Annotation, Transcript};
use std::io::{self, Write};

pub fn write_gff3<W: Write>(
    out: &mut W,
    annotation: &Annotation,
    version: &str,
    command_line: &str,
) -> io::Result<()> {
    writeln!(out, "##gff-version 3")?;
    writeln!(out, "# gffread v{version}")?;
    writeln!(out, "# {command_line}")?;
    for transcript in &annotation.transcripts {
        write_transcript(out, transcript)?;
    }
    Ok(())
}

fn write_transcript<W: Write>(out: &mut W, transcript: &Transcript) -> io::Result<()> {
    write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t.\t{}\t.\tID={}",
        transcript.seqid,
        transcript.source,
        transcript.feature,
        transcript.start,
        transcript.end,
        transcript.strand,
        transcript.id
    )?;
    if let Some(gene_id) = &transcript.gene_id {
        write!(out, ";geneID={gene_id}")?;
    }
    if let Some(gene_name) = &transcript.gene_name {
        write!(out, ";gene_name={gene_name}")?;
    }
    writeln!(out)?;

    for exon in &transcript.exons {
        writeln!(
            out,
            "{}\t{}\texon\t{}\t{}\t.\t{}\t.\tParent={}",
            transcript.seqid, transcript.source, exon.start, exon.end, transcript.strand, transcript.id
        )?;
    }
    for cds in &transcript.cds {
        writeln!(
            out,
            "{}\t{}\tCDS\t{}\t{}\t.\t{}\t{}\tParent={}",
            transcript.seqid, transcript.source, cds.start, cds.end, transcript.strand, cds.phase, transcript.id
        )?;
    }
    Ok(())
}
```

Update `crates/gffread-core/src/lib.rs`:

```rust
pub mod compat;
pub mod emit;
pub mod loader;
pub mod model;
pub mod options;

pub const VERSION: &str = "0.12.9";
```

Replace the `Run` arm in `crates/gffread-cli/src/process.rs` with code that loads the annotation and writes GFF3 when `MainOutput::Gff3` is selected:

```rust
        Ok(CommandMode::Run(options)) => match run_outputs(options) {
            Ok(()) => 0,
            Err(err) => {
                let _ = io::stderr().write_all(err.message.as_bytes());
                err.exit_code
            }
        },
```

Add this helper to `crates/gffread-cli/src/process.rs`:

```rust
use gffread_core::compat::CompatError;
use gffread_core::emit::gff3;
use gffread_core::loader::gff::load_annotation;
use gffread_core::options::{MainOutput, RuntimeOptions};
use std::fs::File;

fn run_outputs(options: RuntimeOptions) -> Result<(), CompatError> {
    let annotation = load_annotation(&options.input)?;
    match options.main_output {
        MainOutput::Gff3 => {
            let output = options
                .output
                .ok_or_else(|| CompatError::new("Error: output file is required in phase-one GFF3 path\n", 1))?;
            let mut file = File::create(&output)
                .map_err(|_| CompatError::new(format!("Error creating file: {}\n", output.display()), 1))?;
            let command_line = format!("gffread {}", options.original_args.join(" "));
            gff3::write_gff3(&mut file, &annotation, gffread_core::VERSION, &command_line)
                .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
            Ok(())
        }
        MainOutput::Gtf | MainOutput::Table => Err(CompatError::new(
            "Error: selected output format is not implemented in this phase-one step\n",
            1,
        )),
    }
}
```

- [ ] **Step 5: Run and close GFF3 byte differences**

Run: `cargo test -p gffread-rs expose_to_simplified_gff3_matches_oracle -- --exact`

Expected: PASS after adjusting transcript feature normalization, attribute selection, and command-line header formatting to match the local C++ oracle. Keep differences in `emit::gff3` or `loader::gff`; do not weaken the oracle assertion.

- [ ] **Step 6: Commit**

```bash
git add crates/gffread-core crates/gffread-cli
git commit -m "feat: emit phase one simplified gff3"
```

## Task 5: Add GTF and Table Output

**Files:**

- Create: `crates/gffread-core/src/emit/gtf.rs`
- Create: `crates/gffread-core/src/emit/table.rs`
- Modify: `crates/gffread-core/src/emit/mod.rs`
- Modify: `crates/gffread-cli/src/process.rs`
- Modify: `crates/gffread-cli/tests/compat_examples.rs`

- [ ] **Step 1: Write failing oracle tests for `-T` and `--table`**

Append to `crates/gffread-cli/tests/compat_examples.rs`:

```rust
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
```

- [ ] **Step 2: Run the format tests and verify they fail**

Run: `cargo test -p gffread-rs --test compat_examples`

Expected: FAIL because `MainOutput::Gtf` and `MainOutput::Table` still return the explicit unsupported-format error.

- [ ] **Step 3: Implement GTF emitter**

Create `crates/gffread-core/src/emit/gtf.rs`:

```rust
use crate::model::{Annotation, Segment, Transcript};
use std::io::{self, Write};

pub fn write_gtf<W: Write>(out: &mut W, annotation: &Annotation) -> io::Result<()> {
    for transcript in &annotation.transcripts {
        write_transcript(out, transcript)?;
    }
    Ok(())
}

fn write_transcript<W: Write>(out: &mut W, transcript: &Transcript) -> io::Result<()> {
    writeln!(
        out,
        "{}\t{}\ttranscript\t{}\t{}\t.\t{}\t.\t{}",
        transcript.seqid,
        transcript.source,
        transcript.start,
        transcript.end,
        transcript.strand,
        attrs_for_transcript(transcript, false)
    )?;
    for exon in &transcript.exons {
        write_segment(out, transcript, exon, "exon")?;
    }
    for cds in &transcript.cds {
        write_segment(out, transcript, cds, "CDS")?;
    }
    Ok(())
}

fn write_segment<W: Write>(
    out: &mut W,
    transcript: &Transcript,
    segment: &Segment,
    feature: &str,
) -> io::Result<()> {
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}",
        transcript.seqid,
        transcript.source,
        feature,
        segment.start,
        segment.end,
        transcript.strand,
        segment.phase,
        attrs_for_transcript(transcript, true)
    )
}

fn attrs_for_transcript(transcript: &Transcript, terminal_semicolon: bool) -> String {
    let gene_id = transcript.gene_id.as_deref().unwrap_or("");
    let gene_name = transcript.gene_name.as_deref().unwrap_or("");
    if terminal_semicolon {
        format!(
            "transcript_id \"{}\"; gene_id \"{}\"; gene_name \"{}\";",
            transcript.id, gene_id, gene_name
        )
    } else {
        format!(
            "transcript_id \"{}\"; gene_id \"{}\"; gene_name \"{}\"",
            transcript.id, gene_id, gene_name
        )
    }
}
```

- [ ] **Step 4: Implement table format parser and table emitter**

Create `crates/gffread-core/src/emit/table.rs`:

```rust
use crate::model::{Annotation, Transcript};
use std::io::{self, Write};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TableField {
    Id,
    Chr,
    Start,
    End,
    Strand,
    Exons,
    Attr(String),
}

pub fn parse_format(format: &str) -> Vec<TableField> {
    format
        .split(|ch| matches!(ch, ',' | ';' | ':' | ' '))
        .filter(|part| !part.is_empty())
        .map(|part| match part {
            "@id" | "ID" | "transcript_id" => TableField::Id,
            "@chr" => TableField::Chr,
            "@start" => TableField::Start,
            "@end" => TableField::End,
            "@strand" => TableField::Strand,
            "@exons" => TableField::Exons,
            other => TableField::Attr(other.to_owned()),
        })
        .collect()
}

pub fn write_table<W: Write>(
    out: &mut W,
    annotation: &Annotation,
    format: &str,
) -> io::Result<()> {
    let fields = parse_format(format);
    for transcript in &annotation.transcripts {
        let row = fields
            .iter()
            .map(|field| value_for(transcript, field))
            .collect::<Vec<_>>()
            .join("\t");
        writeln!(out, "{row}")?;
    }
    Ok(())
}

fn value_for(transcript: &Transcript, field: &TableField) -> String {
    match field {
        TableField::Id => transcript.id.clone(),
        TableField::Chr => transcript.seqid.clone(),
        TableField::Start => transcript.start.to_string(),
        TableField::End => transcript.end.to_string(),
        TableField::Strand => transcript.strand.to_string(),
        TableField::Exons => transcript.exon_list(),
        TableField::Attr(name) => transcript.attrs.get(name).cloned().unwrap_or_default(),
    }
}
```

Update `crates/gffread-core/src/emit/mod.rs`:

```rust
pub mod gff3;
pub mod gtf;
pub mod table;
```

- [ ] **Step 5: Wire GTF and table output in the CLI process layer**

Extend imports in `crates/gffread-cli/src/process.rs`:

```rust
use gffread_core::emit::{gff3, gtf, table};
```

Replace the `match options.main_output` block in `run_outputs` with:

```rust
    match options.main_output {
        MainOutput::Gff3 => {
            let output = options
                .output
                .ok_or_else(|| CompatError::new("Error: output file is required in phase-one GFF3 path\n", 1))?;
            let mut file = File::create(&output)
                .map_err(|_| CompatError::new(format!("Error creating file: {}\n", output.display()), 1))?;
            let command_line = format!("gffread {}", options.original_args.join(" "));
            gff3::write_gff3(&mut file, &annotation, gffread_core::VERSION, &command_line)
                .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
        }
        MainOutput::Gtf => {
            let output = options
                .output
                .ok_or_else(|| CompatError::new("Error: output file is required in phase-one GTF path\n", 1))?;
            let mut file = File::create(&output)
                .map_err(|_| CompatError::new(format!("Error creating file: {}\n", output.display()), 1))?;
            gtf::write_gtf(&mut file, &annotation)
                .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
        }
        MainOutput::Table => {
            let output = options
                .output
                .ok_or_else(|| CompatError::new("Error: output file is required in phase-one table path\n", 1))?;
            let format = options
                .table_format
                .as_deref()
                .ok_or_else(|| CompatError::new("Error: --table requires a format\n", 1))?;
            let mut file = File::create(&output)
                .map_err(|_| CompatError::new(format!("Error creating file: {}\n", output.display()), 1))?;
            table::write_table(&mut file, &annotation, format)
                .map_err(|err| CompatError::new(format!("Error writing output: {err}\n"), 1))?;
        }
    }
    Ok(())
```

- [ ] **Step 6: Run and close GTF/table byte differences**

Run: `cargo test -p gffread-rs --test compat_examples`

Expected: PASS after adjusting transcript feature names, GTF semicolon placement, and attribute aliases to match the oracle. Keep table field parsing in `emit::table::parse_format`.

- [ ] **Step 7: Commit**

```bash
git add crates/gffread-core crates/gffread-cli
git commit -m "feat: emit phase one gtf and table outputs"
```

## Task 6: Add Transcript and CDS FASTA Output

**Files:**

- Create: `crates/gffread-core/src/fasta.rs`
- Modify: `crates/gffread-core/src/lib.rs`
- Modify: `crates/gffread-cli/src/process.rs`
- Modify: `crates/gffread-cli/tests/compat_examples.rs`

- [ ] **Step 1: Write failing oracle tests for `-w` and `-W -x`**

Append to `crates/gffread-cli/tests/compat_examples.rs`:

```rust
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
        .args(["-W", "-x", "transcripts_CDS.fa", "-g", "genome.fa", "annotation.gff"])
        .expected_files(["transcripts_CDS.fa"])
        .assert_matches_oracle(candidate())
        .expect("CDS FASTA must match oracle");
}
```

- [ ] **Step 2: Run FASTA tests and verify they fail**

Run: `cargo test -p gffread-rs --test compat_examples`

Expected: FAIL because no FASTA files are generated.

- [ ] **Step 3: Implement genome loading, splicing, reverse complement, and FASTA wrapping**

Create `crates/gffread-core/src/fasta.rs`:

```rust
use crate::compat::CompatError;
use crate::model::{Annotation, Segment, Transcript};
use std::collections::BTreeMap;
use std::fs;
use std::io::{self, Write};
use std::path::Path;

pub type Genome = BTreeMap<String, Vec<u8>>;

pub fn load_genome(path: &Path) -> Result<Genome, CompatError> {
    let text = fs::read_to_string(path)
        .map_err(|_| CompatError::new(format!("Error: couldn't open genomic sequences file {}\n", path.display()), 1))?;
    let mut genome = Genome::new();
    let mut current_id = String::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            current_id = rest.split_whitespace().next().unwrap_or("").to_owned();
            genome.entry(current_id.clone()).or_default();
        } else if !current_id.is_empty() {
            genome
                .entry(current_id.clone())
                .or_default()
                .extend(line.trim().as_bytes().iter().map(|base| base.to_ascii_uppercase()));
        }
    }
    Ok(genome)
}

pub fn write_transcript_fasta<W: Write>(
    out: &mut W,
    annotation: &Annotation,
    genome: &Genome,
) -> Result<(), CompatError> {
    for transcript in &annotation.transcripts {
        let seq = spliced_sequence(transcript, &transcript.exons, genome)?;
        let defline = transcript_defline(transcript, seq.len());
        write_fasta_record(out, &defline, &seq, false)
            .map_err(|err| CompatError::new(format!("Error writing FASTA: {err}\n"), 1))?;
    }
    Ok(())
}

pub fn write_cds_fasta<W: Write>(
    out: &mut W,
    annotation: &Annotation,
    genome: &Genome,
    write_segments: bool,
) -> Result<(), CompatError> {
    for transcript in &annotation.transcripts {
        if transcript.cds.is_empty() {
            continue;
        }
        let seq = spliced_sequence(transcript, &transcript.cds, genome)?;
        let defline = if write_segments {
            projected_defline(transcript, &transcript.cds)
        } else {
            transcript.id.clone()
        };
        write_fasta_record(out, &defline, &seq, false)
            .map_err(|err| CompatError::new(format!("Error writing FASTA: {err}\n"), 1))?;
    }
    Ok(())
}

pub fn spliced_sequence(
    transcript: &Transcript,
    segments: &[Segment],
    genome: &Genome,
) -> Result<Vec<u8>, CompatError> {
    let chrom = genome.get(&transcript.seqid).ok_or_else(|| {
        CompatError::new(format!("Error: couldn't find genomic sequence {}\n", transcript.seqid), 1)
    })?;
    let mut seq = Vec::new();
    for segment in segments {
        let start = segment.start.saturating_sub(1) as usize;
        let end = segment.end as usize;
        if end > chrom.len() || start >= end {
            return Err(CompatError::new("Error: genomic segment outside sequence bounds\n", 1));
        }
        seq.extend_from_slice(&chrom[start..end]);
    }
    if transcript.strand == '-' {
        reverse_complement(&mut seq);
    }
    Ok(seq)
}

pub fn reverse_complement(seq: &mut Vec<u8>) {
    seq.reverse();
    for base in seq {
        *base = match *base {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            other => other.to_ascii_uppercase(),
        };
    }
}

fn transcript_defline(transcript: &Transcript, seq_len: usize) -> String {
    if transcript.cds.is_empty() {
        transcript.id.clone()
    } else {
        format!("{} CDS=2-{}", transcript.id, seq_len)
    }
}

fn projected_defline(transcript: &Transcript, segments: &[Segment]) -> String {
    let mut offset = 1u64;
    let parts = segments
        .iter()
        .map(|segment| {
            let len = segment.end - segment.start + 1;
            let part = format!("{offset}-{}", offset + len - 1);
            offset += len;
            part
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{} loc:{}({}){}-{} segs:{}",
        transcript.id, transcript.seqid, transcript.strand, transcript.start, transcript.end, parts
    )
}

pub fn write_fasta_record<W: Write>(
    out: &mut W,
    defline: &str,
    seq: &[u8],
    use_star_stop: bool,
) -> io::Result<()> {
    writeln!(out, ">{defline}")?;
    for chunk in seq.chunks(70) {
        if use_star_stop {
            let converted = chunk
                .iter()
                .map(|base| if *base == b'.' { b'*' } else { *base })
                .collect::<Vec<_>>();
            out.write_all(&converted)?;
        } else {
            out.write_all(chunk)?;
        }
        out.write_all(b"\n")?;
    }
    Ok(())
}
```

Update `crates/gffread-core/src/lib.rs`:

```rust
pub mod compat;
pub mod emit;
pub mod fasta;
pub mod loader;
pub mod model;
pub mod options;

pub const VERSION: &str = "0.12.9";
```

- [ ] **Step 4: Wire transcript and CDS FASTA writers**

Add imports to `crates/gffread-cli/src/process.rs`:

```rust
use gffread_core::fasta::{load_genome, write_cds_fasta, write_transcript_fasta};
```

After the main-output `match` in `run_outputs`, add:

```rust
    if options.fasta_outputs.transcript.is_some() || options.fasta_outputs.cds.is_some() {
        let genome_path = options
            .genome
            .as_ref()
            .ok_or_else(|| CompatError::new("Error: -g option is required for options -w/x/y/u/V/N/M !\n", 1))?;
        let genome = load_genome(genome_path)?;

        if let Some(path) = &options.fasta_outputs.transcript {
            let mut file = File::create(path)
                .map_err(|_| CompatError::new(format!("Error creating file: {}\n", path.display()), 1))?;
            write_transcript_fasta(&mut file, &annotation, &genome)?;
        }
        if let Some(path) = &options.fasta_outputs.cds {
            let mut file = File::create(path)
                .map_err(|_| CompatError::new(format!("Error creating file: {}\n", path.display()), 1))?;
            write_cds_fasta(
                &mut file,
                &annotation,
                &genome,
                options.fasta_outputs.write_exon_segments,
            )?;
        }
    }
```

- [ ] **Step 5: Run and close FASTA byte differences**

Run: `cargo test -p gffread-rs --test compat_examples`

Expected: PASS after adjusting CDS coordinate projection, defline formatting, FASTA line width, transcript CDS ranges, and reverse-complement behavior against oracle diffs. Keep sequence logic in `fasta.rs`.

- [ ] **Step 6: Commit**

```bash
git add crates/gffread-core crates/gffread-cli
git commit -m "feat: emit phase one transcript and cds fasta"
```

## Task 7: Add Protein Translation and Combined FASTA Outputs

**Files:**

- Modify: `crates/gffread-core/src/fasta.rs`
- Modify: `crates/gffread-cli/src/process.rs`
- Modify: `crates/gffread-cli/tests/compat_examples.rs`

- [ ] **Step 1: Write failing oracle tests for `-y` and combined `-w -y`**

Append to `crates/gffread-cli/tests/compat_examples.rs`:

```rust
#[test]
fn protein_fasta_matches_oracle() {
    CompatCase::new("example_protein_fasta")
        .in_examples()
        .args(["-y", "transcripts_prot.fa", "-g", "genome.fa", "annotation.gff"])
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
        .expect("combined FASTA outputs must match oracle");
}
```

- [ ] **Step 2: Run protein tests and verify they fail**

Run: `cargo test -p gffread-rs --test compat_examples`

Expected: FAIL because protein FASTA files are not generated.

- [ ] **Step 3: Implement codon translation and protein FASTA writer**

Append to `crates/gffread-core/src/fasta.rs`:

```rust
pub fn write_protein_fasta<W: Write>(
    out: &mut W,
    annotation: &Annotation,
    genome: &Genome,
) -> Result<(), CompatError> {
    for transcript in &annotation.transcripts {
        if transcript.cds.is_empty() {
            continue;
        }
        let cds = spliced_sequence(transcript, &transcript.cds, genome)?;
        let protein = translate(&cds);
        write_fasta_record(out, &transcript.id, protein.as_bytes(), false)
            .map_err(|err| CompatError::new(format!("Error writing FASTA: {err}\n"), 1))?;
    }
    Ok(())
}

pub fn translate(cds: &[u8]) -> String {
    cds.chunks(3)
        .filter(|codon| codon.len() == 3)
        .map(translate_codon)
        .collect()
}

fn translate_codon(codon: &[u8]) -> char {
    match upper_codon(codon).as_str() {
        "TTT" | "TTC" => 'F',
        "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => 'L',
        "ATT" | "ATC" | "ATA" => 'I',
        "ATG" => 'M',
        "GTT" | "GTC" | "GTA" | "GTG" => 'V',
        "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => 'S',
        "CCT" | "CCC" | "CCA" | "CCG" => 'P',
        "ACT" | "ACC" | "ACA" | "ACG" => 'T',
        "GCT" | "GCC" | "GCA" | "GCG" => 'A',
        "TAT" | "TAC" => 'Y',
        "TAA" | "TAG" | "TGA" => '.',
        "CAT" | "CAC" => 'H',
        "CAA" | "CAG" => 'Q',
        "AAT" | "AAC" => 'N',
        "AAA" | "AAG" => 'K',
        "GAT" | "GAC" => 'D',
        "GAA" | "GAG" => 'E',
        "TGT" | "TGC" => 'C',
        "TGG" => 'W',
        "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => 'R',
        "GGT" | "GGC" | "GGA" | "GGG" => 'G',
        _ => 'X',
    }
}

fn upper_codon(codon: &[u8]) -> String {
    codon
        .iter()
        .map(|base| (*base as char).to_ascii_uppercase())
        .collect()
}
```

- [ ] **Step 4: Wire protein output and shared genome loading**

Add `write_protein_fasta` to the FASTA import in `crates/gffread-cli/src/process.rs`:

```rust
use gffread_core::fasta::{
    load_genome, write_cds_fasta, write_protein_fasta, write_transcript_fasta,
};
```

Change the FASTA output condition in `run_outputs` to include protein output:

```rust
    if options.fasta_outputs.transcript.is_some()
        || options.fasta_outputs.cds.is_some()
        || options.fasta_outputs.protein.is_some()
    {
```

Add this block after the CDS writer block:

```rust
        if let Some(path) = &options.fasta_outputs.protein {
            let mut file = File::create(path)
                .map_err(|_| CompatError::new(format!("Error creating file: {}\n", path.display()), 1))?;
            write_protein_fasta(&mut file, &annotation, &genome)?;
        }
```

- [ ] **Step 5: Run and close protein byte differences**

Run: `cargo test -p gffread-rs --test compat_examples`

Expected: PASS after adjusting start phase handling, terminal stop rendering, incomplete trailing codons, and non-coding transcript skipping to match oracle diffs.

- [ ] **Step 6: Commit**

```bash
git add crates/gffread-core/src/fasta.rs crates/gffread-cli
git commit -m "feat: emit phase one protein fasta"
```

## Task 8: Full Phase-One Verification, Make Target, and Function Map

**Files:**

- Create: `docs/compat/function-map.md`
- Modify: `Makefile`
- Modify: `crates/gffread-cli/tests/compat_examples.rs`

- [ ] **Step 1: Add a full-suite test list for the seven example commands**

Ensure `crates/gffread-cli/tests/compat_examples.rs` contains these seven tests and expected files:

```rust
// Required phase-one example cases:
// 1. expose_to_simplified_gff3_matches_oracle -> ann_simple.gff
// 2. gtf_conversion_matches_oracle -> annotation.gtf
// 3. transcript_fasta_matches_oracle -> transcripts.fa
// 4. cds_fasta_with_projected_segments_matches_oracle -> transcripts_CDS.fa
// 5. protein_fasta_matches_oracle -> transcripts_prot.fa
// 6. combined_transcript_and_protein_fasta_match_oracle -> transcripts.fa, transcripts_prot.fa
// 7. table_conversion_matches_oracle -> annotation.tbl
```

- [ ] **Step 2: Run the full oracle suite and verify it passes**

Run: `cargo test -p gffread-rs`

Expected: PASS for version, CLI surface tests, and all seven example tests.

- [ ] **Step 3: Add a Makefile target for Rust compatibility tests**

Append to `Makefile`:

```make

.PHONY : rust-test
rust-test:
	cargo test -p gffread-rs
```

- [ ] **Step 4: Add the function mapping document**

Create `docs/compat/function-map.md`:

```markdown
# GffRead C++ to Rust Function Map

This document tracks phase-one traceability from upstream C++ functions to Rust modules. The C++ source remains the oracle for byte-level compatibility in the current local repository version.

## Process and CLI

- `gffread.cpp::main` maps to `crates/gffread-cli/src/parse.rs` and `crates/gffread-cli/src/process.rs`.
- `gffread.cpp::printGff3Header` maps to `crates/gffread-core/src/emit/gff3.rs::write_gff3`.
- `gffread.cpp::openfw` maps to output-file creation in `crates/gffread-cli/src/process.rs`.

## Parsing and Model

- GCLib `GffLoader::load` behavior needed by phase-one examples maps to `crates/gffread-core/src/loader/gff.rs`.
- GCLib `GffObj` fields needed by phase-one examples map to `crates/gffread-core/src/model.rs::Transcript` and `Segment`.

## Output

- `gff_utils.cpp::printTableData` maps to `crates/gffread-core/src/emit/table.rs::write_table`.
- `gffread.cpp::setTableFormat` maps to `crates/gffread-core/src/emit/table.rs::parse_format`.
- GCLib `GffObj::printGxf` behavior needed for GFF3/GTF examples maps to `crates/gffread-core/src/emit/gff3.rs` and `crates/gffread-core/src/emit/gtf.rs`.
- `gff_utils.cpp::printFasta` maps to `crates/gffread-core/src/fasta.rs::write_fasta_record`.

## Sequence Logic

- GCLib FASTA loading behavior needed by phase-one examples maps to `crates/gffread-core/src/fasta.rs::load_genome`.
- Spliced transcript and CDS extraction maps to `crates/gffread-core/src/fasta.rs::spliced_sequence`.
- Codon translation maps to `crates/gffread-core/src/fasta.rs::translate`.

## Phase-One Boundary

- `gff_utils.cpp::adjust_stopcodon`, `collectLocusData`, `redundantTranscripts`, and advanced filter paths are reserved for later compatibility phases.
- Any CLI option accepted by phase-one parsing but not covered by oracle tests must remain outside the documented complete surface.
```

- [ ] **Step 5: Verify Makefile target**

Run: `make rust-test`

Expected: PASS and equivalent to `cargo test -p gffread-rs`.

- [ ] **Step 6: Commit**

```bash
git add Makefile docs/compat/function-map.md crates/gffread-cli/tests/compat_examples.rs
git commit -m "test: add phase one rust compatibility target"
```

## Final Verification

- [ ] **Step 1: Run formatting**

Run: `cargo fmt --all -- --check`

Expected: PASS. If it fails, run `cargo fmt --all`, inspect the formatting diff, and rerun the check.

- [ ] **Step 2: Run Rust tests**

Run: `cargo test --workspace`

Expected: PASS.

- [ ] **Step 3: Run Rust compatibility target**

Run: `make rust-test`

Expected: PASS.

- [ ] **Step 4: Check git state**

Run: `git status --short`

Expected: no unstaged or untracked files except intentionally ignored build outputs.
