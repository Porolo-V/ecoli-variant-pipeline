# E. coli variant calling pipeline

🇬🇧 English | [🇷🇺 Русский](README.ru.md)

Reproducible short-read variant calling workflow: raw Illumina reads → filtered, annotated variants.
Built on *E. coli* B REL606 (haploid, 4.6 Mb) as a first, deliberately simple case.

## Workflow

```mermaid
flowchart LR
  A[FASTQ] --> B[fastp<br/>trimming]
  B --> C[BWA-MEM2<br/>alignment]
  C --> D[samtools<br/>markdup]
  D --> E[bcftools<br/>mpileup + call]
  E --> F[bcftools filter<br/>QUAL, DP, MQ]
  F --> G[snpEff<br/>annotation]
  G --> H[variants table]
```

Two implementations of the same workflow:
- `scripts/run_all.sh`: plain bash version (also runs FastQC/MultiQC and downloads the data).
- `main.nf` + `nextflow.config`: Nextflow version (6 processes) with a `test` profile.

## Quick start

```bash
mamba env create -f environment.yml
conda activate varcall

# fast test on a 300 kb subset (a few minutes)
nextflow run main.nf -profile test

# full run (downloads data, runs all steps)
bash scripts/run_all.sh
```

## Data

- Reads: SRR2584863 (paired-end Illumina, E. coli evolution experiment), downloaded from ENA.
- Reference: NC_012967.1 (E. coli B str. REL606), FASTA + GenBank from NCBI.

## Results (full sample, bash version)

| Step | Result |
|---|---|
| Reads before / after fastp | 3.10 M / 2.57 M (82.8% kept; Q30 bases 94.4%) |
| Loss reason | 17.2% reads too short after tail trimming, 0.04% low quality |
| Alignment | 100% mapped, 98.2% properly paired |
| Duplicates removed | ~1% |
| Mean depth / breadth | 73x / 99.85% of the genome |
| Variants | 30 raw → 26 (QUAL>30, DP>10) → **22** (MQ>40) |

Of the 22 variants: 13 missense, 3 frameshift, 1 in-frame insertion, 5 up/downstream.
Full table: [`docs/variants_table.tsv`](docs/variants_table.tsv).

<!-- TODO: add docs/multiqc_quality.png and docs/fastp_before_after.png -->

## Filtering decisions

Four variants (positions 263235, 377000, 2407766, 2618472) passed the QUAL/DP filters but had
MQ 22-30, 23-39% of reads with MQ=0 and a mixed ref/alt read split (DP4), which is
unexpected in a haploid organism. This points to repeated regions where reads from several
copies pile up, so a mapping-quality filter (MQ>40) removes them.
The test dataset checks this: position 263235 must be present in the raw VCF and absent after filtering.

## Limitations

- One sample, no orthogonal validation (e.g. IGV inspection, Sanger sequencing).
- Thresholds were chosen by inspecting this sample, not benchmarked.
- Position 2103887 (tandem repeat insertion, depth 28) is kept but has low confidence.
- Indels in homopolymers and repeats are less reliable with short reads.
- The variant table takes the first snpEff `ANN` entry, so intergenic variants can appear as `upstream_gene_variant`.
- `--ploidy 1` is hard-coded in the calling step; other organisms need parameters (planned).
- The Nextflow version was tested on the 300 kb subset only.

## Problems encountered

- Slow connections to GitHub/NCBI/ENA: solved with multi-connection download (`aria2c -x 8`) from ENA.
- FastQC 0.13 uses `--outdir` instead of `-o`.
- The `.gitignore` rule `data/` also ignored `tests/data`; fixed by anchoring it (`/data/`).

## Repository layout

main.nf, nextflow.config Nextflow pipeline and parameters

scripts/run_all.sh bash version of the full workflow

tests/data/ 300 kb subset for the test profile

docs/ variant table, figures

environment.yml conda environment


## Roadmap

Containers (BioContainers), GitHub Actions CI, ploidy and annotation parameters for other organisms.
