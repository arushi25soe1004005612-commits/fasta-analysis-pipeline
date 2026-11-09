#!/usr/bin/env python3
"""FASTA -> stats -> ORFs -> proteins -> GC windows -> report.

Run:  python run_pipeline.py                      (uses the parameters below)
      python run_pipeline.py --input other.fasta  (override input / outdir)
Each step reads the previous step's file from OUTDIR, so every intermediate is kept.
Standard library only.
"""
import argparse, itertools, os, sys

# ============================ PARAMETERS ============================
INPUT      = "data/phix174.fasta"   # input FASTA (first record is analysed)
OUTDIR     = "results"              # all intermediate and final files go here
MIN_ORF_AA = 50                     # minimum ORF length (amino acids, stop excluded)
START      = "ATG"                  # start codon
WINDOW     = 500                    # GC window size (bp)
STEP       = 250                    # GC window step (bp)
# ====================================================================

BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[i] for i, (a, b, c) in enumerate(itertools.product(BASES, BASES, BASES))}
COMP = str.maketrans("ACGTN", "TGCAN")

def revcomp(s): return s.translate(COMP)[::-1]
def translate(s): return "".join(CODON.get(s[i:i+3], "X") for i in range(0, len(s) - 2, 3))
def path(name): return os.path.join(OUTDIR, name)
def wrap(s, n=60): return "\n".join(s[i:i+n] for i in range(0, len(s), n))

def read_fasta(p):
    name, seq = None, []
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            if name is not None: break          # first record only
            name = line[1:].split()[0]
        elif line and name is not None:
            seq.append(line.upper())
    if name is None: sys.exit(f"{p}: no FASTA header found")
    return name, "".join(seq)

# Step 1: validate and basic stats
def step1_stats():
    name, seq = read_fasta(INPUT)
    with open(path("01_clean.fasta"), "w") as f: f.write(f">{name}\n{wrap(seq)}\n")
    gc = 100 * (seq.count("G") + seq.count("C")) / max(1, sum(seq.count(b) for b in "ACGT"))
    with open(path("01_stats.tsv"), "w") as f:
        f.write("id\tlength\tGC_percent\tA\tC\tG\tT\tother\n")
        f.write(f"{name}\t{len(seq)}\t{gc:.2f}\t" + "\t".join(str(seq.count(b)) for b in "ACGT")
                + f"\t{len(seq) - sum(seq.count(b) for b in 'ACGT')}\n")

# Step 2: find ORFs in all six frames (start codon to stop, forward and reverse strand)
def step2_orfs():
    name, seq = read_fasta(path("01_clean.fasta"))
    L, rows = len(seq), []
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        for frame in range(3):
            i = frame
            while i + 3 <= L:
                if s[i:i+3] == START:
                    j = i
                    while j + 3 <= L and CODON.get(s[j:j+3]) != "*": j += 3
                    if j + 3 <= L and (j - i) // 3 >= MIN_ORF_AA:       # stop codon found
                        a, b = (i + 1, j + 3) if strand == "+" else (L - j - 2, L - i)
                        rows.append((strand, frame + 1, a, b, (j - i) // 3, s[i:j+3]))
                    i = j + 3
                else:
                    i += 3
    rows.sort(key=lambda r: -r[4])
    with open(path("02_orfs.tsv"), "w") as f:
        f.write("orf_id\tstrand\tframe\tstart\tend\tlength_aa\n")
        for n, r in enumerate(rows, 1): f.write(f"ORF{n:02d}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\n")
    with open(path("02_orfs.fasta"), "w") as f:
        for n, r in enumerate(rows, 1): f.write(f">ORF{n:02d} {r[0]} {r[2]}-{r[3]}\n{wrap(r[5])}\n")

# Step 3: translate the ORFs
def step3_proteins():
    recs, hdr = [], None
    for line in open(path("02_orfs.fasta")):
        line = line.strip()
        if line.startswith(">"): hdr = line[1:]; recs.append([hdr, ""])
        elif line: recs[-1][1] += line
    with open(path("03_proteins.faa"), "w") as f:
        for h, s in recs: f.write(f">{h}\n{wrap(translate(s).rstrip('*'))}\n")

# Step 4: GC content in sliding windows
def step4_gc_windows():
    _, seq = read_fasta(path("01_clean.fasta"))
    with open(path("04_gc_windows.tsv"), "w") as f:
        f.write("start\tend\tGC_percent\n")
        for a in range(0, len(seq) - WINDOW + 1, STEP):
            w = seq[a:a + WINDOW]
            f.write(f"{a + 1}\t{a + WINDOW}\t{100 * (w.count('G') + w.count('C')) / WINDOW:.2f}\n")

# Step 5: report built only from the earlier files
def step5_report():
    stats = open(path("01_stats.tsv")).read().splitlines()[1].split("\t")
    orfs = [l.split("\t") for l in open(path("02_orfs.tsv")).read().splitlines()[1:]]
    gcw = [(l.split("\t")) for l in open(path("04_gc_windows.tsv")).read().splitlines()[1:]]
    lo = min(gcw, key=lambda r: float(r[2])); hi = max(gcw, key=lambda r: float(r[2]))
    with open(path("05_report.md"), "w") as f:
        f.write(f"# Pipeline report\n\nInput: `{INPUT}`\n\n")
        f.write(f"Parameters: MIN_ORF_AA={MIN_ORF_AA}, START={START}, WINDOW={WINDOW}, STEP={STEP}\n\n")
        f.write(f"- Sequence `{stats[0]}`: {stats[1]} bp, GC {stats[2]}%\n")
        f.write(f"- ORFs of at least {MIN_ORF_AA} aa: {len(orfs)} "
                f"({sum(o[1] == '+' for o in orfs)} forward, {sum(o[1] == '-' for o in orfs)} reverse)\n")
        f.write(f"- Lowest-GC window: {lo[0]}-{lo[1]} ({lo[2]}%); highest: {hi[0]}-{hi[1]} ({hi[2]}%)\n\n")
        f.write("| ORF | strand | frame | start | end | aa |\n|---|---|---|---|---|---|\n")
        for o in orfs[:10]: f.write("| " + " | ".join(o) + " |\n")

STEPS = [step1_stats, step2_orfs, step3_proteins, step4_gc_windows, step5_report]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--input"); ap.add_argument("--outdir")
    a = ap.parse_args()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))   # paths are relative to the repo
    if a.input: INPUT = a.input
    if a.outdir: OUTDIR = a.outdir
    os.makedirs(OUTDIR, exist_ok=True)
    for s in STEPS:
        print(f"running {s.__name__} ...")
        s()
    print(f"done. Files in {OUTDIR}/:", ", ".join(sorted(os.listdir(OUTDIR))))
