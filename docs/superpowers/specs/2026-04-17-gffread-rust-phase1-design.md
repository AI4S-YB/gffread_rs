# GffRead Rust Phase 1 Design

## Goal

Rewrite the current repository version of `gffread` in Rust with an oracle-first process. Phase 1 establishes a Rust executable, a byte-level compatibility harness, and enough implementation to match the current sample commands and CLI surface cases against the C++ `gffread v0.12.9` source in this workspace.

## Scope

The oracle is the C++ `gffread` built from the current repository contents. Phase 1 targets the current local Linux/WSL environment only. The compatibility contract includes exit code, stdout, stderr, help/version text, error text, and all generated output files.

Phase 1 includes:

- Create a Rust workspace and `gffread-rs` executable entry point.
- Reproduce process-level behavior for `-h`, `--help`, `--version`, invalid option combinations, missing required options, missing files, stdout, stderr, and exit codes.
- Build an oracle-first harness that runs the C++ binary and Rust binary with the same inputs and compares captured artifacts byte-for-byte.
- Match all commands represented by `gffread/run_tests.sh`, including generated file bytes.
- Add module skeletons and mapping documentation so later phases can trace core Rust behavior back to upstream C++ functions.

Phase 1 explicitly does not include:

- Full parity for every public CLI option.
- Full implementation of clustering, advanced filtering, streaming, all BED/TLF/GTF edge cases, complex comment retention, or every sequence-validation path unless needed by the phase-one cases.
- Cross-platform byte parity.
- Treating a parsed or accepted option as complete unless it is covered by oracle tests and documented as implemented.

## Recommended Approach

Use an oracle-first compatibility workflow.

The C++ implementation remains the source of truth. Each compatibility test first runs the C++ `gffread`, captures exit code, stdout, stderr, and generated files, then runs the Rust binary and compares every captured artifact. This approach is slower to set up than a direct rewrite, but it gives a hard correctness signal for the stated byte-level parity goal.

Rejected alternatives:

- Direct whole-program translation before tests: initially fast, but leaves late-stage byte diffs hard to localize.
- Mixed implementation and harness construction: produces earlier Rust output, but weakens the feedback loop and increases rework.

## Architecture

The Rust implementation will be organized as a workspace with clear compatibility boundaries.

### `crates/gffread-cli`

Owns process-level behavior:

- Argument parsing.
- Help and version text.
- CLI validation ordering.
- File open checks.
- Conversion from internal errors to C++-compatible stderr text and exit codes.

This layer should intentionally follow the order of checks in `gffread.cpp::main`, because byte-level behavior depends on which condition is reported first.

### `crates/gffread-core`

Owns domain behavior and is split by responsibility:

- `options`: normalized runtime options after CLI parsing.
- `loader`: GFF/GTF/BED/TLF input parsing and comment/header handling.
- `model`: genes, transcripts, exons, CDS segments, loci, and related attributes.
- `filters`: ID, range, junction, CDS, pseudo, length, intron, and coding/non-coding filters.
- `cluster`: merge, cluster-only, redundancy, locus collection, and coverage behavior.
- `fasta`: transcript/CDS extraction, translation, stop codon adjustment, and FASTA formatting.
- `emit`: GFF3, GTF, BED, TLF, table, and junction output.
- `compat`: byte-sensitive compatibility details such as formatting, ordering, warnings, and upstream quirks.

Core algorithms should preserve traceability to upstream functions, while Rust-specific I/O, error propagation, ownership, and resource management can be idiomatic.

### `crates/gffread-test-harness`

Owns compatibility testing:

- Build or locate the C++ oracle binary from the local `gffread` source tree.
- Run the C++ oracle and Rust candidate in isolated temporary directories.
- Capture exit code, stdout, stderr, and generated output files.
- Compare artifacts byte-for-byte.
- Emit useful diffs for failures.

Text diffs should use unified diff output. Non-text or difficult byte diffs should report the first differing byte offset and a compact hex context.

### `fixtures`

Phase 1 should reuse `gffread/examples` as the primary fixture source. Additional minimized fixtures should be added only when they isolate a specific compatibility behavior that the sample commands do not cover.

### `docs/compat`

Documents upstream-to-Rust traceability:

- Which C++ function maps to which Rust module/function.
- Which mappings are exact algorithm ports.
- Which mappings intentionally use Rust-specific structure while preserving observable behavior.
- Which options or branches are present but not yet phase-one complete.

## Data Flow

The phase-one runtime flow is:

`argv -> CLI parsing -> runtime options -> input reader -> parsed records -> transcript model -> output writers -> captured artifacts -> oracle comparison`

All output writers should write to `Write` implementations or byte buffers and avoid platform-dependent line handling. Output formatting is part of the compatibility surface, not a presentation detail.

## Function Mapping Policy

Use a hybrid function-to-function policy:

- Core algorithm paths should be traceable to the C++ implementation.
- Peripheral concerns such as ownership, file handles, structured errors, and test orchestration can be idiomatic Rust.
- Behavior that exists only for C++ compatibility should live in clearly named compatibility modules rather than being hidden inside general-purpose APIs.

Initial mapping targets:

- `gffread.cpp::main` -> `gffread-cli` parsing and validation flow.
- `loadIDlist` -> `ids::load_id_list`.
- `loadSeqInfo` -> `seqinfo::load`.
- `loadRefTable` -> `refs::load_rename_table`.
- `getAttrList` -> `attrs::parse_filter_list`.
- `setTableFormat` -> `emit::table::parse_format`.
- `printGff3Header` -> `emit::gff3::write_header`.
- `printGSeqHeader` -> `emit::gff3::write_sequence_region`.
- `processGffComment` -> `loader::comments::process_gff_comment`.
- `printGffObj` -> `emit::gff::write_object`.
- `printAsTable` and `printTableData` -> `emit::table::write_record`.
- `printFasta` -> `fasta::write_record`.
- `adjust_stopcodon` -> `cds::adjust_stop_codon`.
- `collectLocusData` -> `cluster::collect_locus_data`.
- `redundantTranscripts` -> `cluster::redundant_transcripts`.

## Testing Strategy

Testing has three layers.

### Process Compatibility Tests

Run the same command through C++ and Rust and compare:

- Exit code.
- Stdout bytes.
- Stderr bytes.
- Generated file names.
- Generated file bytes.

The first required command set is the seven commands encoded in `gffread/run_tests.sh`.

### CLI Surface Tests

Add oracle-backed cases for:

- `-h`.
- `--help`.
- `--version`.
- Missing `-g` when sequence output or validation requires it.
- Invalid option combinations such as `--sort-by` with `--sort-alpha`.
- Missing input file or output creation failure.

### Minimized Function-Path Tests

Add small fixtures as needed for:

- Table field parsing.
- ID list parsing.
- Sequence info parsing.
- Reference rename table parsing.
- GFF3 header generation.
- FASTA wrapping and defline formatting.
- Stop codon adjustment boundaries.

These tests should stay oracle-backed when observable output is involved.

## Error Handling

Rust internals should use structured errors. Public process output must go through a compatibility renderer that controls:

- Exact message text.
- Newline behavior.
- Whether text goes to stdout or stderr.
- Exit code.
- Ordering of validation failures.

No user-visible error should be emitted by directly formatting an arbitrary Rust error if byte parity is expected.

## Completion Criteria

Phase 1 is complete when:

- The Rust workspace builds.
- The `gffread-rs` executable can be launched by the harness.
- Oracle-backed tests for the sample commands and CLI surface cases pass with matching exit code, stdout, stderr, and generated file bytes.
- Any remaining unsupported options are explicitly documented and not presented as complete.
- `docs/compat` contains enough C++-to-Rust mapping information to continue later phases without redesigning the project.

