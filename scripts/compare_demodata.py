#!/usr/bin/env python3
import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


CORE_FEATURES = {"gene", "mRNA", "transcript", "exon", "CDS"}


def run(cmd, cwd):
    return subprocess.run(
        [str(part) for part in cmd],
        cwd=cwd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )


def run_to_files(cmd, cwd, stdout_path, stderr_path):
    with stdout_path.open("wb") as stdout_handle, stderr_path.open("wb") as stderr_handle:
        completed = subprocess.run(
            [str(part) for part in cmd],
            cwd=cwd,
            stdout=stdout_handle,
            stderr=stderr_handle,
            check=False,
        )
    return completed.returncode


def clean_id(value):
    return value.split(",")[0].strip()


def parse_attrs(raw):
    attrs = {}
    for part in raw.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
        else:
            continue
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def discover_cases(demodata):
    cases = []
    for gff in sorted(demodata.glob("*/*.gff3")):
        stem = gff.name
        if ".representtaive." in stem:
            continue
        fasta_candidates = list(gff.parent.glob("*.genome.fasta")) + list(
            gff.parent.glob("*.genome.fna")
        )
        if not fasta_candidates:
            continue
        cases.append((gff.parent.name, gff, fasta_candidates[0]))
    return cases


def pick_transcripts(gff, limit):
    selected = []
    genes = set()
    with gff.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            if fields[2] not in {"mRNA", "transcript"}:
                continue
            attrs = parse_attrs(fields[8])
            transcript_id = attrs.get("ID") or attrs.get("transcript_id")
            if not transcript_id:
                continue
            selected.append(clean_id(transcript_id))
            parent = attrs.get("Parent") or attrs.get("geneID") or attrs.get("gene_id")
            if parent:
                genes.add(clean_id(parent))
            if len(selected) >= limit:
                break
    return set(selected), genes


def extract_subset(src, dest, limit):
    if limit <= 0:
        link_or_copy(src, dest)
        return

    transcripts, genes = pick_transcripts(src, limit)
    if not transcripts:
        raise RuntimeError(f"no transcript records found in {src}")

    kept = 0
    with src.open() as input_handle, dest.open("w") as output_handle:
        for line in input_handle:
            if line.startswith("#"):
                output_handle.write(line)
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in CORE_FEATURES:
                continue
            attrs = parse_attrs(fields[8])
            keep = False
            if fields[2] == "gene":
                feature_id = attrs.get("ID") or attrs.get("gene_id")
                keep = feature_id is not None and clean_id(feature_id) in genes
            elif fields[2] in {"mRNA", "transcript"}:
                feature_id = attrs.get("ID") or attrs.get("transcript_id")
                keep = feature_id is not None and clean_id(feature_id) in transcripts
            else:
                parent = attrs.get("Parent") or attrs.get("transcript_id")
                keep = parent is not None and clean_id(parent) in transcripts
            if keep:
                output_handle.write(line)
                kept += 1
    if kept == 0:
        raise RuntimeError(f"subset extraction produced no feature lines for {src}")


def compare_bytes(path_a, path_b):
    if not path_a.exists() or not path_b.exists():
        return False
    return subprocess.run(
        ["cmp", "-s", str(path_a), str(path_b)],
        check=False,
    ).returncode == 0


def first_diff(path_a, path_b):
    if not path_a.exists() and not path_b.exists():
        return "both files are missing"
    if not path_a.exists():
        return f"oracle file missing, rust size {path_b.stat().st_size}"
    if not path_b.exists():
        return f"rust file missing, oracle size {path_a.stat().st_size}"
    a = path_a.read_bytes() if path_a.exists() else b""
    b = path_b.read_bytes() if path_b.exists() else b""
    if a == b:
        return "identical"
    limit = min(len(a), len(b))
    idx = next((i for i in range(limit) if a[i] != b[i]), limit)
    a_line = a[:idx].count(b"\n") + 1
    b_line = b[:idx].count(b"\n") + 1
    return f"first byte diff at {idx}, oracle line {a_line}, rust line {b_line}, sizes {len(a)} vs {len(b)}"


def sanitize(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def link_or_copy(src, dest):
    try:
        os.link(src, dest)
    except OSError:
        try:
            os.symlink(src, dest)
        except OSError:
            shutil.copy2(src, dest)


def run_case(case_name, gff, fasta, workdir, oracle, rust, sample_size):
    case_dir = workdir / sanitize(case_name)
    shutil.rmtree(case_dir, ignore_errors=True)
    oracle_dir = case_dir / "oracle"
    rust_dir = case_dir / "rust"
    oracle_dir.mkdir(parents=True)
    rust_dir.mkdir(parents=True)

    subset = case_dir / "sample.gff3"
    extract_subset(gff, subset, sample_size)
    fasta_arg = str(fasta.resolve())

    if any(suffix in fasta.name for suffix in (".fasta", ".fna", ".fa")):
        run([rust, "-w", os.devnull, "-g", fasta_arg, str(subset.resolve())], workdir)

    tasks = [
        ("gff3", ["sample.gff3"], ["stdout"]),
        ("gtf", ["-T", "sample.gff3"], ["stdout"]),
        ("exon_fa", ["-w", "exons.fa", "-g", fasta_arg, "sample.gff3"], ["exons.fa"]),
        ("cds_fa", ["-x", "cds.fa", "-g", fasta_arg, "sample.gff3"], ["cds.fa"]),
        ("protein_fa", ["-y", "protein.fa", "-g", fasta_arg, "sample.gff3"], ["protein.fa"]),
        (
            "combined_fa",
            [
                "-w",
                "combined_exons.fa",
                "-x",
                "combined_cds.fa",
                "-y",
                "combined_protein.fa",
                "-g",
                fasta_arg,
                "sample.gff3",
            ],
            ["combined_exons.fa", "combined_cds.fa", "combined_protein.fa"],
        ),
    ]

    failures = []
    for task_name, args, expected_files in tasks:
        task_oracle = oracle_dir / task_name
        task_rust = rust_dir / task_name
        task_oracle.mkdir()
        task_rust.mkdir()
        link_or_copy(subset, task_oracle / "sample.gff3")
        link_or_copy(subset, task_rust / "sample.gff3")

        oracle_stdout = task_oracle / "stdout"
        oracle_stderr = task_oracle / "stderr"
        rust_stdout = task_rust / "stdout"
        rust_stderr = task_rust / "stderr"

        oracle_code = run_to_files([oracle, *args], task_oracle, oracle_stdout, oracle_stderr)
        rust_code = run_to_files([rust, *args], task_rust, rust_stdout, rust_stderr)

        if oracle_code != rust_code:
            failures.append(
                f"{case_name}/{task_name}: exit {oracle_code} != {rust_code}"
            )
        if not compare_bytes(oracle_stdout, rust_stdout):
            failures.append(
                f"{case_name}/{task_name}: stdout differs ({first_diff(oracle_stdout, rust_stdout)})"
            )
        if not compare_bytes(oracle_stderr, rust_stderr):
            failures.append(
                f"{case_name}/{task_name}: stderr differs ({first_diff(oracle_stderr, rust_stderr)})"
            )
        for expected in expected_files:
            if expected == "stdout":
                continue
            oracle_file = task_oracle / expected
            rust_file = task_rust / expected
            if not compare_bytes(oracle_file, rust_file):
                failures.append(
                    f"{case_name}/{task_name}: {expected} differs ({first_diff(oracle_file, rust_file)})"
                )

        if compare_bytes(oracle_stdout, rust_stdout):
            oracle_stdout.unlink(missing_ok=True)
            rust_stdout.unlink(missing_ok=True)
        if compare_bytes(oracle_stderr, rust_stderr):
            oracle_stderr.unlink(missing_ok=True)
            rust_stderr.unlink(missing_ok=True)

    return failures


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--demodata", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, default=Path("target/demodata-compat"))
    parser.add_argument("--oracle", type=Path, default=Path("./gffread"))
    parser.add_argument("--rust", type=Path, default=Path("target/debug/gffread-rs"))
    parser.add_argument("--sample-size", type=int, default=10)
    parser.add_argument("--case", action="append", default=[])
    args = parser.parse_args()

    cases = discover_cases(args.demodata)
    if args.case:
        wanted = set(args.case)
        cases = [case for case in cases if case[0] in wanted]

    if not cases:
        print("no cases discovered", file=sys.stderr)
        return 2

    args.oracle = args.oracle.resolve()
    args.rust = args.rust.resolve()
    args.workdir.mkdir(parents=True, exist_ok=True)
    failures = []
    for case_name, gff, fasta in cases:
        print(f"CASE {case_name}: {gff.name} / {fasta.name}", flush=True)
        failures.extend(
            run_case(
                case_name,
                gff,
                fasta,
                args.workdir,
                args.oracle,
                args.rust,
                args.sample_size,
            )
        )

    if failures:
        print("\nFAILURES:")
        for failure in failures:
            print(f"- {failure}")
        return 1

    print(f"all {len(cases)} cases matched")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
