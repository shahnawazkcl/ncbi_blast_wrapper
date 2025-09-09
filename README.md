# NCBI BLASTp Wrapper (R) — WIP

![Status](https://img.shields.io/badge/status-work--in--progress-yellow)
![GitHub last commit](https://img.shields.io/github/last-commit/shahnawazkcl/ncbi_blast_wrapper)
[![MIT License](https://img.shields.io/badge/License-MIT%20License-orange)](https://opensource.org/license/mit)
![Commit](https://img.shields.io/github/commit-activity/t/shahnawazkcl/ncbi_blast_wrapper)


Small R pipeline that runs offline **BLASTp** pairwise alignments for curated pairs of protein sequences.

> **Project status:** early setup / scaffold. Expect breaking changes.

## Requirements
- **NCBI BLAST+** (tested with 2.12.0+)  
- **R** ≥ 4.0 with package: `Biostrings`  
- (Optional) RStudio

_The original scaffold targeted BLAST+ 2.12.0+, R 4.0.3._  
See `library(Biostrings)` in the script for package needs.

## Repo layout

.  
├─ `ncbi\_blastp\_wrap.R`       # main pipeline  
├─ `prepare\_input.sh`          # builds ./data from \*.fasta in repo root  
├─ `transcript\_type\_info.csv`  # mapping of pairs: name, prin, alt  
└─ data/                     # example multi-FASTA files  


## Input expectations

- `transcript_type_info.csv` with columns:

|name|prin|alt
-----|-----|---
ESRRB\_1|ENST00000512784\_domains.fasta|ENST00000505752\_domains.fasta
ESRRB\_2|ENST00000512784\_domains.fasta|ENST00000644823\_domains.fasta

- `./data/<name>.fasta` multi-FASTA files whose sequence headers contain the transcript IDs used above (e.g., `>zf-C4_1_ENST00000512784`).

## Quick start

1) Install BLAST+ and R deps (`Biostrings`).
2) Prepare data (either):
 - Put your multi-FASTA files under `./data`, **or**
 - Place `*.fasta` in repo root and run:
   ```bash
   bash prepare_input.sh
   ```
3) Configure `ncbi_blastp_wrap.R` if needed:
 - `path_to_domain_files <- "data/"` (default)
4) Run:
 ```bash
 Rscript ncbi_blastp_wrap.R
````

## What it does (pipeline)

* Splits each multi-FASTA into per-sequence files and tags them as **query** or **subject** based on `transcript_type_info.csv`.
* For each `name`, finds the matching query/subject pair and runs:

  ```bash
  blastp -query <query.fasta> -subject <subject.fasta> -outfmt 0
  ```
* Saves text outputs under `./alignment/` (folder name may appear as `alingment/` in early versions).

## Outputs

* For each domain/family `name`, a folder is created containing:

  * `alignment/Alignment_<query>_<subject>_out.txt` (pairwise BLASTp report, `-outfmt 0`).

## Roadmap / known issues (early stage)

* Define or replace `open_input_files()`; fix `save_metafile()` variable names/scope.
* Normalize output folder to `alignment/`.
* Add argument parsing (input dir, CSV path, outdir, `-outfmt`).
* Add unit tests and CI for R/BLAST+ presence.
* Example notebook / vignette.

## License

TBD.
