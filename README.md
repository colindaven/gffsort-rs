# gffsort-rs

A Rust implementation of a GFF3 sorter designed for tabix-compatible ordering.

This project is a reimplementation of the Perl script `gff3sort.pl` from the original `tao-bioinfo/gff3sort` project:

https://github.com/tao-bioinfo/gff3sort/blob/main/gff3sort.pl

It keeps the same core behavior for sorting GFF3 records by chromosome and start position while preserving correct ordering for features that share the same chromosome and coordinate, especially in `--precise` mode where parent features are emitted before children.

## Why this exists

Tabix requires GFF3 files to be sorted by chromosome and start position. A naive sort can place parent features after their children when multiple records share the same coordinate, which breaks downstream display in tools such as JBrowse.

This tool follows the same intent as the Perl reference script: group records by chromosome and start position, then keep feature ordering stable or topologically sorted when needed.

## Features

- Sorts GFF3 records by chromosome and start position
- Supports chromosome ordering modes:
  - `alphabet` (default is `natural`)
  - `natural`
  - `original`
- Supports `--precise` mode to order parent features before children when needed
- Preserves pragma lines such as `##sequence-region` and other `#...` header metadata
- Ignores inline FASTA blocks unless `--extract-fasta` is used
- Optionally extracts inline FASTA to a separate `.fasta` file
- Supports writing output to a file with `-o/--output`

## Usage

Build:

```bash
cargo build --release
```

Run to stdout:

```bash
cargo run --release -- input.gff3 > sorted.gff3
```

Write directly to a file:

```bash
cargo run --release -- -o output.gff input.gff3
```

Or with explicit options:

```bash
cargo run --release -- --chr_order natural --precise -o output.gff input.gff3
```

Extract FASTA alongside sorting:

```bash
cargo run --release -- --extract-fasta -o sorted.gff3 input.gff3
```

Read from standard input:

```bash
cat input.gff3 | cargo run --release -- --chr_order natural --precise -o output.gff -
```

## Command-line options

```text
-o, --output <OUTPUT>
    Write sorted GFF3 output to FILE instead of stdout.

--chr_order <alphabet|natural|original>
    Chromosome sort order. Default: natural

--precise
    Topologically sort records at the same chromosomal position so
    parent features appear before their children.

--extract-fasta
    Extract inline FASTA sequence records to a separate .fasta file.

INPUT
    Optional GFF3 input file path. Use '-' or omit to read from stdin.
```

## Compatibility notes

This Rust version mirrors the behavior and intent of the original Perl implementation, but it is not a byte-for-byte port. It aims to provide the same practical sorting behavior for GFF3 files used in tabix indexing workflows.

## License

This project is distributed under the GPL-3.0 license.
