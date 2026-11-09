# Pipeline report

Input: `data/phix174.fasta`

Parameters: MIN_ORF_AA=50, START=ATG, WINDOW=500, STEP=250

- Sequence `NC_001422.1`: 5386 bp, GC 44.76%
- ORFs of at least 50 aa: 22 (16 forward, 6 reverse)
- Lowest-GC window: 1-500 (41.00%); highest: 501-1000 (48.80%)

| ORF | strand | frame | start | end | aa |
|---|---|---|---|---|---|
| ORF01 | + | 2 | 1001 | 2284 | 427 |
| ORF02 | + | 3 | 2931 | 3917 | 328 |
| ORF03 | + | 1 | 3076 | 3684 | 202 |
| ORF04 | - | 1 | 1094 | 1651 | 185 |
| ORF05 | + | 1 | 2395 | 2922 | 175 |
| ORF06 | + | 3 | 390 | 848 | 152 |
| ORF07 | - | 1 | 1772 | 2167 | 131 |
| ORF08 | + | 1 | 3946 | 4263 | 105 |
| ORF09 | - | 3 | 3921 | 4211 | 96 |
| ORF10 | + | 1 | 568 | 843 | 91 |
