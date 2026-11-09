# FASTA analysis pipeline

One command runs a five-step pipeline on a FASTA file and keeps every intermediate file.

```bash
python run_pipeline.py
```

Requirements: Python 3.7+ only (standard library, nothing to install).
Input: `data/phix174.fasta` (phiX174 complete genome, NCBI RefSeq NC_001422.1, included in the repo).
Output: everything is written to `results/`.

To use another file: `python run_pipeline.py --input my.fasta --outdir my_results`

## Parameters
All parameters are in one block at the top of `run_pipeline.py`:

| Parameter | Default | Meaning |
|---|---|---|
| `INPUT` | `data/phix174.fasta` | input FASTA (first record is analysed) |
| `OUTDIR` | `results` | folder for all output files |
| `MIN_ORF_AA` | 50 | minimum ORF length in amino acids |
| `START` | ATG | start codon |
| `WINDOW` / `STEP` | 500 / 250 | GC sliding-window size and step (bp) |

## Steps and intermediate files
Each step reads the file the previous step wrote, so any step can be inspected on its own.

| Step | Does | Writes |
|---|---|---|
| 1 | validate input, length, GC, base counts | `01_clean.fasta`, `01_stats.tsv` |
| 2 | find ORFs in all six frames (start codon to stop codon) | `02_orfs.tsv`, `02_orfs.fasta` |
| 3 | translate the ORFs | `03_proteins.faa` |
| 4 | GC content in sliding windows | `04_gc_windows.tsv` |
| 5 | summary report | `05_report.md` |

## Results on phiX174 (5,386 bp, GC 44.76%)
22 ORFs of at least 50 aa (16 forward, 6 reverse). The longest ones match known phiX174 genes:
gene F (major capsid protein, 1001-2284, 427 aa), gene H (2931-3917, 328 aa), gene G (2395-2922, 175 aa) and gene E (lysis protein, 568-843, 91 aa).
Full output: `results/05_report.md`.

## Decisions and limits
- **Standard library only**, so a clean checkout runs with no installs.
- **Data is in the repo**, so the run is reproducible without downloading anything. The script works from any folder because paths are relative to the repo.
- **ORF definition:** ATG to the first in-frame stop, both strands. Nested ATGs inside a longer ORF in the same frame are not reported separately.
- **Linear sequence:** phiX174 is circular, so an ORF that runs across the end of the sequence (for example gene A, which starts at 3981) is not found whole. This is a known simplification.
- The ORFs are predictions from sequence only. I compared the longest ones with known gene positions but did not run a gene-annotation tool.
