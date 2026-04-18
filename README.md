## GffRead 

GFF/GTF utility providing format conversions, filtering, FASTA sequence 
extraction and more.

More details and usage examples can be found in the paper [DOI: 10.12688/f1000research.23297.1](http://dx.doi.org/10.12688/f1000research.23297.1) which can be also used to cite this software.

The official webpage with download packages for this utility can be found online here: 
 http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread

Use `gffread -h` to see the command line usage options.

## Installation
Building this program from source requires the [GCLib](../../../gclib) source code 
library. The `make` command should automatically fetch the latest gclib version from the repository if no `../gclib` directory is found.

```
  cd /some/build/dir
  git clone https://github.com/gpertea/gffread
  cd gffread
  make release
```
This should create the **gffread** binary in the current directory.

## Rust rewrite status

This repository also contains an in-progress Rust rewrite of `gffread`.
The Rust command line binary is built from the Cargo workspace as
`gffread-rs`, while the repository-local C++ `./gffread` binary is used as the
compatibility oracle.

Build the optimized Rust binary with:

```bash
cargo build -p gffread-rs --release
```

The release binary is written to:

```text
target/release/gffread-rs
```

## Compatibility evaluation

The current parity target is function-to-function algorithmic compatibility and
byte-to-byte output compatibility against the current repository version of
`./gffread`.

The full demodata validation compares these output classes for every case:

- GFF3 conversion on stdout
- GTF conversion on stdout
- exon/transcript FASTA extraction
- CDS FASTA extraction
- protein FASTA extraction
- combined exon/CDS/protein FASTA extraction

The validated full-input species set is:

- `Amborella_trichopoda`
- `Arabidopsis_thaliana`
- `Citrus_sinensis`
- `Malus_domestica`
- `Oryza_sativa`
- `Solanum_lycopersicum`

The full validation command used for the latest run was:

```bash
python3 scripts/compare_demodata.py \
  --demodata /mnt/j/ClaudeCodeDev/gffread_rs/demodata \
  --workdir /tmp/gffread-demodata-full-final \
  --sample-size 0 \
  --oracle ./gffread \
  --rust target/release/gffread-rs
```

`--sample-size 0` means the full input file for each species is used. The
`/tmp` work directory is intentional for WSL/Linux runs because it avoids
DrvFS-mounted filesystem stalls that can make large FASTA comparisons appear
much slower than the actual binaries.

Latest full-run result:

```text
all 6 cases matched
demodata_full_final elapsed=768.59 user=244.18 sys=33.33 maxrss=1923020
```

This means that, for the demodata evaluation scope above, the Rust binary's
GFF3/GTF conversions and exon/CDS/protein FASTA outputs matched the repository
C++ `./gffread` oracle byte-for-byte.

Additional verification commands used before publishing the Rust parity work:

```bash
cargo fmt --check
python3 -m unittest scripts.tests.test_compare_demodata
cargo test -p gffread-rs --test compat_examples -- --nocapture
cargo test -p gffread-rs --test compat_fasta_options -- --nocapture
```

The latest local compatibility test counts were:

```text
compat_examples: 48 passed
compat_fasta_options: 12 passed
compare_demodata unit tests: 4 passed
```

## Demodata comparison script

The demodata comparison helper is:

```text
scripts/compare_demodata.py
```

It runs each selected case with both the C++ oracle and the Rust binary, writes
stdout/stderr and generated FASTA files under the requested work directory, and
uses byte comparisons for all expected outputs. When a mismatch is found, it
reports the task, stream or file, first differing byte, corresponding line
number, and output sizes.

Example single-species full run:

```bash
python3 scripts/compare_demodata.py \
  --demodata /mnt/j/ClaudeCodeDev/gffread_rs/demodata \
  --workdir /tmp/gffread-demodata-oryza \
  --sample-size 0 \
  --case Oryza_sativa \
  --oracle ./gffread \
  --rust target/release/gffread-rs
```
