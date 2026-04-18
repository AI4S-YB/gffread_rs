use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::ffi::{OsStr, OsString};
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Mutex, OnceLock};

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

    pub fn assert_matches_oracle(&self, candidate: impl AsRef<Path>) -> Result<(), Box<dyn Error>> {
        let oracle_run = {
            let _guard = oracle_run_lock()
                .lock()
                .expect("oracle test harness lock should not be poisoned");
            let oracle = ensure_oracle_binary()?;
            self.run_command(&oracle)?
        };
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

        if let Err(failure) = compare_output_file_sets(
            &self.name,
            &oracle_run.output_files,
            &candidate_run.output_files,
        ) {
            failures.push(failure);
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
        let files_before = capture_regular_files(&workdir, &workdir)?;

        let output = Command::new(program)
            .args(&self.args)
            .current_dir(&workdir)
            .output()?;
        let files_after = capture_regular_files(&workdir, &workdir)?;

        Ok(RunResult {
            status_code: output.status.code().unwrap_or(-1),
            stdout: output.stdout,
            stderr: output.stderr,
            output_files: changed_regular_file_paths(&files_before, &files_after),
            files: capture_expected_files(&workdir, &self.expected_files)?,
        })
    }
}

fn oracle_run_lock() -> &'static Mutex<()> {
    static ORACLE_RUN_LOCK: OnceLock<Mutex<()>> = OnceLock::new();
    ORACLE_RUN_LOCK.get_or_init(|| Mutex::new(()))
}

struct RunResult {
    status_code: i32,
    stdout: Vec<u8>,
    stderr: Vec<u8>,
    output_files: BTreeSet<PathBuf>,
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
    ensure_oracle_binary_with_make(&repo_root(), Path::new("make"))
}

fn ensure_oracle_binary_with_make(
    repo_root: &Path,
    make: &Path,
) -> Result<PathBuf, Box<dyn Error>> {
    let oracle = repo_root.join("gffread");
    let mut command = Command::new(make);
    command.arg("-C").arg(repo_root);

    if let Some(gcldir) = discover_gclib_dir(repo_root) {
        command.arg(format!("GCLDIR={}", gcldir.display()));
    }

    let status = command.arg("gffread").status()?;
    if !status.success() {
        return Err("failed to build C++ oracle with make".into());
    }
    if !oracle.is_file() {
        return Err("C++ oracle binary missing after make".into());
    }
    Ok(oracle)
}

fn discover_gclib_dir(repo_root: &Path) -> Option<PathBuf> {
    let mut candidates = Vec::new();
    if let Some(parent) = repo_root.parent() {
        candidates.push(parent.join("gclib"));
    }
    candidates.push(repo_root.join(".worktrees").join("gclib"));

    candidates.into_iter().find_map(|candidate| {
        let has_headers = candidate.join("GBase.h").is_file() && candidate.join("gff.h").is_file();
        if has_headers {
            Some(candidate.canonicalize().unwrap_or(candidate))
        } else {
            None
        }
    })
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

fn capture_regular_files(
    root: &Path,
    dir: &Path,
) -> Result<BTreeMap<PathBuf, Vec<u8>>, Box<dyn Error>> {
    let mut captured = BTreeMap::new();
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        let file_type = entry.file_type()?;
        if file_type.is_dir() {
            captured.extend(capture_regular_files(root, &path)?);
        } else if file_type.is_file() {
            captured.insert(path.strip_prefix(root)?.to_path_buf(), fs::read(path)?);
        }
    }
    Ok(captured)
}

fn changed_regular_file_paths(
    before: &BTreeMap<PathBuf, Vec<u8>>,
    after: &BTreeMap<PathBuf, Vec<u8>>,
) -> BTreeSet<PathBuf> {
    after
        .iter()
        .filter_map(|(path, after_bytes)| match before.get(path) {
            Some(before_bytes) if before_bytes == after_bytes => None,
            Some(_) | None => Some(path.clone()),
        })
        .collect()
}

fn compare_output_file_sets(
    case_name: &str,
    oracle_files: &BTreeSet<PathBuf>,
    candidate_files: &BTreeSet<PathBuf>,
) -> Result<(), String> {
    if oracle_files == candidate_files {
        return Ok(());
    }

    let mut message = format!("{case_name}: output file set differed");
    let _ = write!(&mut message, "\noracle: {}", format_path_set(oracle_files));
    let _ = write!(
        &mut message,
        "\ncandidate: {}",
        format_path_set(candidate_files)
    );
    Err(message)
}

fn format_path_set(paths: &BTreeSet<PathBuf>) -> String {
    if paths.is_empty() {
        return "[]".to_string();
    }

    let joined = paths
        .iter()
        .map(|path| path.display().to_string())
        .collect::<Vec<_>>()
        .join(", ");
    format!("[{joined}]")
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

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(unix)]
    use std::os::unix::fs::PermissionsExt;

    #[cfg(not(windows))]
    const REBUILD_MOCK_MAKE: &str = "#!/bin/sh\nprintf invoked > \"$2/gffread\"\n";
    #[cfg(windows)]
    const REBUILD_MOCK_MAKE: &str =
        "@echo off\r\n> \"%~2\\gffread\" <nul set /p \"=invoked\"\r\nexit /b 0\r\n";

    #[cfg(not(windows))]
    const GCLDIR_MOCK_MAKE: &str = "#!/bin/sh\nfor arg in \"$@\"; do\n  case \"$arg\" in\n    GCLDIR=*) printf '%s' \"${arg#GCLDIR=}\" > \"$2/gcldir.txt\" ;;\n  esac\ndone\nprintf invoked > \"$2/gffread\"\n";
    #[cfg(windows)]
    const GCLDIR_MOCK_MAKE: &str = "@echo off\r\nsetlocal EnableDelayedExpansion\r\nset \"arg=%~3\"\r\nif /I \"!arg:~0,7!\"==\"GCLDIR=\" set \"gcldir=!arg:~7!\"\r\nif defined gcldir > \"%~2\\gcldir.txt\" <nul set /p \"=!gcldir!\"\r\n> \"%~2\\gffread\" <nul set /p \"=invoked\"\r\nexit /b 0\r\n";

    fn write_mock_make(bin_dir: &Path, script: &str) -> PathBuf {
        #[cfg(windows)]
        let make_path = bin_dir.join("make.cmd");
        #[cfg(not(windows))]
        let make_path = bin_dir.join("make");

        fs::write(&make_path, script).expect("make script should be written");

        #[cfg(unix)]
        {
            let mut permissions = fs::metadata(&make_path)
                .expect("make script metadata should exist")
                .permissions();
            permissions.set_mode(0o755);
            fs::set_permissions(&make_path, permissions).expect("make script should be executable");
        }

        make_path
    }

    #[test]
    fn different_output_filenames_with_same_count_are_reported() {
        let oracle_files = BTreeSet::from([PathBuf::from("oracle.txt")]);
        let candidate_files = BTreeSet::from([PathBuf::from("candidate.txt")]);

        let failure = compare_output_file_sets("case", &oracle_files, &candidate_files)
            .expect_err("different output filenames must fail");

        assert!(failure.contains("output file set differed"));
        assert!(failure.contains("oracle.txt"));
        assert!(failure.contains("candidate.txt"));
    }

    #[test]
    fn run_command_reports_existing_file_modified_in_place() {
        let case = CompatCase::new("modifies existing fixture")
            .in_examples()
            .args(["-c", "printf modified > README.md"]);

        let run = case
            .run_command(Path::new("sh"))
            .expect("shell command should run");

        assert!(
            run.output_files.contains(Path::new("README.md")),
            "modified existing fixture should be reported as changed output; got {}",
            format_path_set(&run.output_files)
        );
    }

    #[test]
    fn ensure_oracle_binary_rebuilds_even_when_binary_exists() {
        let tempdir = TempDir::new().expect("tempdir should be created");
        let repo_root = tempdir.path().join("repo");
        fs::create_dir(&repo_root).expect("repo root should be created");

        fs::write(repo_root.join("gffread"), b"stale").expect("stale binary should be created");

        let bin_dir = tempdir.path().join("bin");
        fs::create_dir(&bin_dir).expect("bin dir should be created");
        let make_path = write_mock_make(&bin_dir, REBUILD_MOCK_MAKE);

        let oracle = ensure_oracle_binary_with_make(&repo_root, &make_path)
            .expect("oracle binary should be ensured");

        assert_eq!(oracle, repo_root.join("gffread"));
        assert_eq!(
            fs::read_to_string(oracle).expect("oracle should be readable"),
            "invoked",
            "build helper should refresh the binary even when it already exists"
        );
    }

    #[test]
    fn ensure_oracle_binary_passes_repo_worktrees_gcldir_for_primary_checkout() {
        let tempdir = TempDir::new().expect("tempdir should be created");
        let repo_root = tempdir.path().join("repo");
        fs::create_dir(&repo_root).expect("repo root should be created");

        let gclib_dir = repo_root.join(".worktrees").join("gclib");
        fs::create_dir_all(&gclib_dir).expect("gclib dir should be created");
        fs::write(gclib_dir.join("GBase.h"), "").expect("GBase header should exist");
        fs::write(gclib_dir.join("gff.h"), "").expect("gff header should exist");

        let bin_dir = tempdir.path().join("bin");
        fs::create_dir(&bin_dir).expect("bin dir should be created");
        let make_path = write_mock_make(&bin_dir, GCLDIR_MOCK_MAKE);

        ensure_oracle_binary_with_make(&repo_root, &make_path)
            .expect("oracle binary should be ensured");

        assert_eq!(
            fs::read_to_string(repo_root.join("gcldir.txt"))
                .expect("make invocation should record the discovered gclib path"),
            gclib_dir
                .canonicalize()
                .expect("gclib path should resolve")
                .display()
                .to_string()
        );
    }
}
