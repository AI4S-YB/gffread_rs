import importlib.util
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "compare_demodata.py"


def load_compare_demodata():
    spec = importlib.util.spec_from_file_location("compare_demodata", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def write_fake_gffread(path: Path):
    path.write_text(
        textwrap.dedent(
            """\
            #!/usr/bin/env python3
            import pathlib
            import sys

            args = sys.argv[1:]
            input_path = pathlib.Path(args[-1])
            text = input_path.read_text()

            if "-T" in args or all(flag not in args for flag in ("-w", "-x", "-y")):
                sys.stdout.write(text)

            sys.stderr.write("note\\n")

            for index, arg in enumerate(args):
                if arg in {"-w", "-x", "-y"}:
                    pathlib.Path(args[index + 1]).write_text(f"{arg}:{input_path.name}\\n")
            """
        )
    )
    path.chmod(0o755)


def write_fixture_gff(path: Path):
    path.write_text(
        textwrap.dedent(
            """\
            ##gff-version 3
            chr1\tsrc\tgene\t1\t12\t.\t+\t.\tID=gene1
            chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=tx1;Parent=gene1
            chr1\tsrc\texon\t1\t12\t.\t+\t.\tParent=tx1
            chr1\tsrc\tCDS\t1\t12\t.\t+\t0\tParent=tx1
            """
        )
    )


def write_fixture_fasta(path: Path):
    path.write_text(">chr1\nACGTACGTACGT\n")


class CompareDemodataTests(unittest.TestCase):
    def setUp(self):
        self.module = load_compare_demodata()
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.gff = self.root / "case.gff3"
        self.fasta = self.root / "case.genome.fa"
        self.fake = self.root / "fake_gffread.py"
        write_fixture_gff(self.gff)
        write_fixture_fasta(self.fasta)
        write_fake_gffread(self.fake)

    def tearDown(self):
        self.tempdir.cleanup()

    def test_run_to_files_streams_process_output_to_paths(self):
        stdout_path = self.root / "stdout.txt"
        stderr_path = self.root / "stderr.txt"

        returncode = self.module.run_to_files(
            [self.fake, str(self.gff)],
            self.root,
            stdout_path,
            stderr_path,
        )

        self.assertEqual(returncode, 0)
        self.assertEqual(stdout_path.read_text(), self.gff.read_text())
        self.assertEqual(stderr_path.read_text(), "note\n")

    def test_run_discards_process_output_for_prewarm(self):
        wrapper = self.root / "wrapper.py"
        wrapper.write_text(
            textwrap.dedent(
                f"""\
                import importlib.util
                from pathlib import Path

                script_path = Path({str(SCRIPT_PATH)!r})
                spec = importlib.util.spec_from_file_location("compare_demodata", script_path)
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                result = module.run([Path({str(self.fake)!r}), Path({str(self.gff)!r})], Path({str(self.root)!r}))
                print(result.returncode)
                """
            )
        )

        completed = subprocess.run(
            [sys.executable, str(wrapper)],
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(completed.returncode, 0)
        self.assertEqual(completed.stdout, "0\n")
        self.assertEqual(completed.stderr, "")

    def test_run_case_accepts_zero_sample_size_as_full_input(self):
        workdir = self.root / "work"

        failures = self.module.run_case(
            "case",
            self.gff,
            self.fasta,
            workdir,
            self.fake,
            self.fake,
            0,
        )

        self.assertEqual(failures, [])
        shared_input = workdir / "case" / "sample.gff3"
        self.assertTrue(shared_input.exists())
        self.assertTrue(os.path.samefile(shared_input, self.gff))

    def test_task_inputs_reuse_shared_sample_file(self):
        workdir = self.root / "work"

        failures = self.module.run_case(
            "case",
            self.gff,
            self.fasta,
            workdir,
            self.fake,
            self.fake,
            1,
        )

        self.assertEqual(failures, [])
        shared_input = workdir / "case" / "sample.gff3"
        self.assertTrue(shared_input.exists())
        self.assertTrue(
            os.path.samefile(workdir / "case" / "oracle" / "gff3" / "sample.gff3", shared_input)
        )
        self.assertTrue(
            os.path.samefile(
                workdir / "case" / "rust" / "combined_fa" / "sample.gff3", shared_input
            )
        )


if __name__ == "__main__":
    unittest.main()
