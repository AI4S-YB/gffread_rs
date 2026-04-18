# STAR-Style CI and Release Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the single Rust workflow with STAR_rs-style CI and release automation for `gffread-rs`.

**Architecture:** Split automation into `.github/workflows/ci.yml` for validation and `.github/workflows/release.yml` for tag/manual releases. Preserve the existing gffread oracle compatibility test split in CI, while using the STAR_rs release pattern for multi-target binary artifacts.

**Tech Stack:** GitHub Actions, Rust stable toolchain, `dtolnay/rust-toolchain`, `Swatinem/rust-cache`, `taiki-e/upload-rust-binary-action`, GitHub CLI.

---

### Task 1: Replace Single Rust Workflow With CI Workflow

**Files:**
- Delete: `.github/workflows/rust.yml`
- Create: `.github/workflows/ci.yml`

- [ ] **Step 1: Delete the old workflow**

Remove `.github/workflows/rust.yml` so the new CI workflow is the only validation workflow for branch pushes and pull requests.

- [ ] **Step 2: Create `.github/workflows/ci.yml`**

Write this exact workflow:

```yaml
name: CI

"on":
  pull_request:
    branches: ["master"]
  push:
    branches: ["master"]

concurrency:
  group: ci-${{ github.ref }}
  cancel-in-progress: true

env:
  CARGO_TERM_COLOR: always

jobs:
  fmt:
    name: Format
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
        with:
          components: rustfmt
      - run: cargo fmt --all -- --check

  clippy:
    name: Clippy
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
        with:
          components: clippy
      - uses: Swatinem/rust-cache@v2
        with:
          key: clippy
      - run: cargo clippy --workspace --all-targets -- -D warnings

  test:
    name: Test (${{ matrix.os }})
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, macos-latest, windows-latest]
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - uses: Swatinem/rust-cache@v2
        with:
          key: test-${{ matrix.os }}
      - name: Build workspace
        run: cargo build --workspace --verbose
      - name: Run Rust-only tests
        run: |
          cargo test -p gffread-core --verbose
          cargo test -p gffread-test-harness --verbose
      - name: Run CLI smoke tests
        run: cargo test -p gffread-rs --test version_smoke --verbose
      - name: Prepare C++ oracle
        if: runner.os != 'Windows'
        run: |
          git clone --depth 1 https://github.com/gpertea/gclib.git .worktrees/gclib
          make GCLDIR="$PWD/.worktrees/gclib" gffread
      - name: Run oracle compatibility tests
        if: runner.os != 'Windows'
        run: |
          cargo test -p gffread-rs --test version_oracle --verbose
          cargo test -p gffread-rs --test compat_cli --verbose
          cargo test -p gffread-rs --test compat_examples --verbose
          cargo test -p gffread-rs --test compat_fasta_options --verbose

  msrv:
    name: MSRV (1.88)
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@1.88.0
      - uses: Swatinem/rust-cache@v2
        with:
          key: msrv
      - run: cargo check --workspace --all-targets
```

- [ ] **Step 3: Validate YAML shape**

Run:

```bash
python3 - <<'PY'
from pathlib import Path
for path in [Path(".github/workflows/ci.yml")]:
    seen = set()
    for i, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip() or line.startswith((" ", "-")):
            continue
        key = line.split(":", 1)[0].strip().strip("\"'")
        if key in seen:
            raise SystemExit(f"{path}: duplicate top-level key {key!r} at line {i}")
        seen.add(key)
    print(f"{path}: no duplicate top-level keys: {sorted(seen)}")
PY
```

Expected: the script prints the four top-level keys without duplicate errors.

- [ ] **Step 4: Commit CI workflow**

```bash
git add .github/workflows/ci.yml .github/workflows/rust.yml
git commit -m "ci: adopt STAR-style validation workflow"
```

### Task 2: Add STAR-Style Release Workflow

**Files:**
- Create: `.github/workflows/release.yml`

- [ ] **Step 1: Create `.github/workflows/release.yml`**

Write this exact workflow:

```yaml
name: Release

"on":
  push:
    tags: ["v*"]
  workflow_dispatch:
    inputs:
      tag:
        description: "Tag to build (e.g. v0.1.0)"
        required: true
        default: "v0.1.0"

permissions:
  contents: write

env:
  CARGO_TERM_COLOR: always

jobs:
  create-release:
    name: Create draft release
    runs-on: ubuntu-latest
    outputs:
      tag: ${{ steps.tag.outputs.tag }}
    steps:
      - uses: actions/checkout@v4
      - name: Resolve tag
        id: tag
        env:
          INPUT_TAG: ${{ inputs.tag }}
        run: |
          if [ "$GITHUB_EVENT_NAME" = "workflow_dispatch" ]; then
            TAG="$INPUT_TAG"
          else
            TAG="${GITHUB_REF#refs/tags/}"
          fi
          echo "tag=$TAG" >> "$GITHUB_OUTPUT"
      - name: Create draft release if missing
        env:
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
          TAG: ${{ steps.tag.outputs.tag }}
        run: |
          if gh release view "$TAG" >/dev/null 2>&1; then
            echo "Release $TAG already exists, skipping create"
          else
            gh release create "$TAG" --draft --generate-notes --title "$TAG" --target "$GITHUB_SHA"
          fi

  build:
    name: Build (${{ matrix.target }})
    needs: create-release
    runs-on: ${{ matrix.os }}
    timeout-minutes: 30
    strategy:
      fail-fast: false
      matrix:
        include:
          - target: x86_64-unknown-linux-gnu
            os: ubuntu-latest
          - target: x86_64-unknown-linux-musl
            os: ubuntu-latest
            apt: musl-tools
          - target: aarch64-apple-darwin
            os: macos-14
          - target: x86_64-pc-windows-msvc
            os: windows-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
        with:
          targets: ${{ matrix.target }}
      - uses: Swatinem/rust-cache@v2
        with:
          key: release-${{ matrix.target }}
      - name: Install apt deps
        if: matrix.apt != ''
        run: |
          sudo apt-get update
          sudo apt-get install -y ${{ matrix.apt }}
      - uses: taiki-e/upload-rust-binary-action@v1
        with:
          bin: gffread-rs
          target: ${{ matrix.target }}
          archive: gffread-rs-${{ needs.create-release.outputs.tag }}-${{ matrix.target }}
          checksum: sha256
          tar: unix
          zip: windows
          ref: refs/tags/${{ needs.create-release.outputs.tag }}
          token: ${{ secrets.GITHUB_TOKEN }}

  publish-release:
    name: Publish release
    needs: [create-release, build]
    runs-on: ubuntu-latest
    steps:
      - name: Publish
        env:
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
          GH_REPO: ${{ github.repository }}
          TAG: ${{ needs.create-release.outputs.tag }}
        run: gh release edit "$TAG" --draft=false
```

- [ ] **Step 2: Validate YAML shape**

Run the duplicate top-level key script from Task 1 against both `.github/workflows/ci.yml` and `.github/workflows/release.yml`.

Expected: both files report no duplicate top-level keys.

- [ ] **Step 3: Commit release workflow**

```bash
git add .github/workflows/release.yml
git commit -m "ci: add STAR-style release workflow"
```

### Task 3: Document CI and Release Usage

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Add README automation section**

Add this section near the end of `README.md`:

```markdown
## CI and releases

GitHub Actions runs the CI workflow on pushes and pull requests targeting `master`. The workflow checks formatting, runs clippy, builds the workspace, runs Rust-only tests on Linux/macOS/Windows, and runs C++ oracle compatibility tests on Linux and macOS.

Releases are built by pushing a `v*` tag, for example:

```bash
git tag v0.1.0
git push ai4s v0.1.0
```

The Release workflow can also be started manually from GitHub Actions with a tag input. Release assets include platform-specific `gffread-rs` binaries for Linux GNU, Linux musl, macOS arm64, and Windows MSVC, plus SHA-256 checksums.
```

- [ ] **Step 2: Verify Markdown fences**

Run:

```bash
python3 - <<'PY'
from pathlib import Path
text = Path("README.md").read_text()
if text.count("```") % 2:
    raise SystemExit("README.md has an unmatched fenced code block")
print("README.md fenced code blocks are balanced")
PY
```

Expected: `README.md fenced code blocks are balanced`.

- [ ] **Step 3: Commit README docs**

```bash
git add README.md
git commit -m "docs: describe CI and release workflow"
```

### Task 4: Full Local Verification and Push

**Files:**
- Verify: `.github/workflows/ci.yml`
- Verify: `.github/workflows/release.yml`
- Verify: `README.md`

- [ ] **Step 1: Run formatting**

```bash
cargo fmt --check
```

Expected: exit code 0.

- [ ] **Step 2: Run clippy**

```bash
cargo clippy --workspace --all-targets -- -D warnings
```

Expected: exit code 0.

- [ ] **Step 3: Run build and core tests**

```bash
cargo build --workspace --verbose
cargo test -p gffread-core --verbose
cargo test -p gffread-test-harness --verbose
cargo test -p gffread-rs --test version_smoke --verbose
```

Expected: all commands exit code 0.

- [ ] **Step 4: Run oracle tests on local Unix environment**

```bash
make GCLDIR="$PWD/.worktrees/gclib" gffread
cargo test -p gffread-rs --test version_oracle --verbose
cargo test -p gffread-rs --test compat_cli --verbose
cargo test -p gffread-rs --test compat_examples --verbose
cargo test -p gffread-rs --test compat_fasta_options --verbose
```

Expected: all commands exit code 0.

- [ ] **Step 5: Push to `ai4s/master`**

```bash
git push ai4s master
```

Expected: remote branch updates successfully.

- [ ] **Step 6: Verify GitHub Actions**

```bash
gh run list --repo AI4S-YB/gffread_rs --limit 5 --json databaseId,headSha,status,conclusion,workflowName,displayTitle,createdAt,url
```

Watch the newest CI run for the pushed commit:

```bash
gh run watch <run-id> --repo AI4S-YB/gffread_rs --exit-status
```

Expected: the CI workflow completes successfully on Linux, macOS, and Windows.
