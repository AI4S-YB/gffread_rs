use std::error::Error;
use std::ffi::{OsStr, OsString};
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use similar::TextDiff;
use tempfile::TempDir;

pub struct CompatCase {
    name: String,
    args: Vec<OsString>,
    in_examples: bool,
    expected_files: Vec<PathBuf>,
}

impl CompatCase {
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            args: Vec::new(),
            in_examples: false,
            expected_files: Vec::new(),
        }
    }

    pub fn args<I, S>(mut self, args: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        self.args = args
            .into_iter()
            .map(|arg| arg.as_ref().to_os_string())
            .collect();
        self
    }

    pub fn in_examples(mut self) -> Self {
        self.in_examples = true;
        self
    }

    pub fn expected_files<I, P>(mut self, paths: I) -> Self
    where
        I: IntoIterator<Item = P>,
        P: AsRef<Path>,
    {
        self.expected_files = paths
            .into_iter()
            .map(|path| path.as_ref().to_path_buf())
            .collect();
        self
    }

    pub fn assert_matches_oracle(
        &self,
        candidate: impl AsRef<Path>,
    ) -> Result<(), Box<dyn Error>> {
        let oracle = ensure_oracle_binary()?;
        let oracle_run = self.run_command(&oracle)?;
        let candidate_run = self.run_command(candidate.as_ref())?;

        let mut failures = Vec::new();

        if oracle_run.status_code != candidate_run.status_code {
            failures.push(format!(
                "{}: exit status differed: oracle={} candidate={}",
                self.name, oracle_run.status_code, candidate_run.status_code
            ));
        }

        if oracle_run.stdout != candidate_run.stdout {
            failures.push(diff_bytes(
                &format!("{} stdout", self.name),
                &oracle_run.stdout,
                &candidate_run.stdout,
            ));
        }

        if oracle_run.stderr != candidate_run.stderr {
            failures.push(diff_bytes(
                &format!("{} stderr", self.name),
                &oracle_run.stderr,
                &candidate_run.stderr,
            ));
        }

        if oracle_run.files.len() != candidate_run.files.len() {
            failures.push(format!(
                "{}: output file count differed: oracle={} candidate={}",
                self.name,
                oracle_run.files.len(),
                candidate_run.files.len()
            ));
        }

        for expected in &self.expected_files {
            let oracle_bytes = oracle_run.files.get(expected);
            let candidate_bytes = candidate_run.files.get(expected);
            match (oracle_bytes, candidate_bytes) {
                (Some(oracle_bytes), Some(candidate_bytes)) if oracle_bytes != candidate_bytes => {
                    failures.push(diff_bytes(
                        &format!("{} file {}", self.name, expected.display()),
                        oracle_bytes,
                        candidate_bytes,
                    ));
                }
                (Some(_), None) => failures.push(format!(
                    "{}: candidate did not create expected file {}",
                    self.name,
                    expected.display()
                )),
                (None, Some(_)) => failures.push(format!(
                    "{}: oracle did not create expected file {}",
                    self.name,
                    expected.display()
                )),
                (Some(_), Some(_)) | (None, None) => {}
            }
        }

        if failures.is_empty() {
            Ok(())
        } else {
            Err(failures.join("\n\n").into())
        }
    }

    fn run_command(&self, program: &Path) -> Result<RunResult, Box<dyn Error>> {
        let tempdir = TempDir::new()?;
        let workdir = if self.in_examples {
            let examples_dir = tempdir.path().join("examples");
            copy_dir_recursive(&repo_root().join("examples"), &examples_dir)?;
            examples_dir
        } else {
            tempdir.path().to_path_buf()
        };

        let output = Command::new(program)
            .args(&self.args)
            .current_dir(&workdir)
            .output()?;

        Ok(RunResult {
            status_code: output.status.code().unwrap_or(-1),
            stdout: output.stdout,
            stderr: output.stderr,
            files: capture_expected_files(&workdir, &self.expected_files)?,
        })
    }
}

struct RunResult {
    status_code: i32,
    stdout: Vec<u8>,
    stderr: Vec<u8>,
    files: Vec<(PathBuf, Vec<u8>)>,
}

trait FileMapExt {
    fn get(&self, path: &Path) -> Option<&Vec<u8>>;
}

impl FileMapExt for Vec<(PathBuf, Vec<u8>)> {
    fn get(&self, path: &Path) -> Option<&Vec<u8>> {
        self.iter()
            .find(|(candidate, _)| candidate == path)
            .map(|(_, bytes)| bytes)
    }
}

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .canonicalize()
        .expect("repo root must resolve")
}

fn ensure_oracle_binary() -> Result<PathBuf, Box<dyn Error>> {
    let oracle = repo_root().join("gffread");
    if oracle.is_file() {
        return Ok(oracle);
    }

    let status = Command::new("make")
        .arg("-C")
        .arg(repo_root())
        .status()?;
    if !status.success() {
        return Err("failed to build C++ oracle with make".into());
    }
    if !oracle.is_file() {
        return Err("C++ oracle binary missing after make".into());
    }
    Ok(oracle)
}

fn copy_dir_recursive(src: &Path, dst: &Path) -> Result<(), Box<dyn Error>> {
    fs::create_dir_all(dst)?;
    for entry in fs::read_dir(src)? {
        let entry = entry?;
        let src_path = entry.path();
        let dst_path = dst.join(entry.file_name());
        let file_type = entry.file_type()?;
        if file_type.is_dir() {
            copy_dir_recursive(&src_path, &dst_path)?;
        } else if file_type.is_file() {
            fs::copy(src_path, dst_path)?;
        }
    }
    Ok(())
}

fn capture_expected_files(
    workdir: &Path,
    expected_files: &[PathBuf],
) -> Result<Vec<(PathBuf, Vec<u8>)>, Box<dyn Error>> {
    let mut captured = Vec::with_capacity(expected_files.len());
    for path in expected_files {
        let full_path = workdir.join(path);
        if full_path.exists() {
            captured.push((path.clone(), fs::read(full_path)?));
        }
    }
    Ok(captured)
}

fn diff_bytes(label: &str, oracle: &[u8], candidate: &[u8]) -> String {
    let oracle_text = String::from_utf8_lossy(oracle);
    let candidate_text = String::from_utf8_lossy(candidate);
    let diff = TextDiff::from_lines(&oracle_text, &candidate_text);
    let mut rendered = String::new();
    let _ = writeln!(&mut rendered, "{} differed:", label);
    let _ = write!(
        &mut rendered,
        "{}",
        diff.unified_diff().header("oracle", "candidate")
    );
    rendered
}
