# STAR-Style CI and Release Design

## Goal

Adopt the CI and release structure used by `AI4S-YB/STAR_rs` while preserving the current `gffread_rs` compatibility test coverage.

## Scope

This change is limited to GitHub Actions automation and minimal documentation. It does not rename the Rust binary, change package metadata, or alter gffread behavior.

## Workflow Structure

Replace the current single `.github/workflows/rust.yml` workflow with two workflows:

- `.github/workflows/ci.yml` for validation on branch pushes and pull requests.
- `.github/workflows/release.yml` for tag-driven and manually dispatched release builds.

The old `rust.yml` workflow will be deleted to avoid duplicate CI runs.

## CI Workflow

The CI workflow mirrors the `STAR_rs` layout:

- Trigger on `push` and `pull_request` for `master`.
- Use concurrency keyed by ref, with in-progress runs cancelled for the same ref.
- Run separate jobs for formatting, clippy, tests, and MSRV checking.
- Use `dtolnay/rust-toolchain@stable` and `Swatinem/rust-cache@v2`.

The `test` matrix runs on `ubuntu-latest`, `macos-latest`, and `windows-latest`. It keeps the existing `gffread_rs` test split:

- All platforms run workspace build, Rust-only tests, and the `version_smoke` CLI smoke test.
- Linux and macOS additionally clone `gpertea/gclib`, build the C++ oracle with `GCLDIR`, and run oracle compatibility tests.
- Windows skips the C++ oracle tests because the current oracle harness and upstream C++ build path are Unix-oriented.

The MSRV job follows the `STAR_rs` pattern with Rust `1.88.0` unless the build proves incompatible; if incompatible, the workflow should use the lowest version that can compile this workspace.

## Release Workflow

The release workflow mirrors the `STAR_rs` release flow:

- Trigger on pushed tags matching `v*`.
- Support `workflow_dispatch` with a required `tag` input.
- Use `permissions.contents: write`.
- Create a draft GitHub release if it does not already exist.
- Build and upload release binaries for:
  - `x86_64-unknown-linux-gnu`
  - `x86_64-unknown-linux-musl`
  - `aarch64-apple-darwin`
  - `x86_64-pc-windows-msvc`
- Publish the draft release only after all build jobs complete.

The release binary remains `gffread-rs`. Archive names follow:

```text
gffread-rs-${tag}-${target}
```

Each uploaded archive includes a SHA-256 checksum.

## Release Implementation Details

Use `taiki-e/upload-rust-binary-action@v1`, matching `STAR_rs`, with:

- `bin: gffread-rs`
- `target: ${{ matrix.target }}`
- `archive: gffread-rs-${{ needs.create-release.outputs.tag }}-${{ matrix.target }}`
- `checksum: sha256`
- `tar: unix`
- `zip: windows`
- `ref: refs/tags/${{ needs.create-release.outputs.tag }}`
- `token: ${{ secrets.GITHUB_TOKEN }}`

Linux musl builds install `musl-tools` before building.

## Documentation

Add a short README section describing:

- CI runs on `master` pushes and pull requests.
- Releases are created by pushing a `v*` tag or manually dispatching the Release workflow.
- Release assets contain platform-specific `gffread-rs` binaries and SHA-256 checksums.

## Validation

Before pushing, run:

```bash
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings
cargo build --workspace --verbose
cargo test -p gffread-core --verbose
cargo test -p gffread-test-harness --verbose
cargo test -p gffread-rs --test version_smoke --verbose
make GCLDIR="$PWD/.worktrees/gclib" gffread
cargo test -p gffread-rs --test version_oracle --verbose
cargo test -p gffread-rs --test compat_cli --verbose
cargo test -p gffread-rs --test compat_examples --verbose
cargo test -p gffread-rs --test compat_fasta_options --verbose
```

After pushing, verify the CI workflow run for the pushed commit completes successfully on Linux, macOS, and Windows.
