# scHLA-demux

`scHLA-demux` is a Snakemake workflow for donor demultiplexing of 10x single-cell data starting from raw FASTQ files. The workflow runs:

1. `cellranger multi` to generate the filtered barcode matrix and sample BAM
2. `souporcell` to infer donor clusters from the gene-expression BAM
3. `arcasHLA` on each souporcell singlet cluster to genotype HLA alleles
4. An optional HLA-matching step to map souporcell clusters back to known donor HLA annotations
5. A publishing step that collects donor-specific outputs into a ready-to-use demultiplexed folder

It is a cleaned-up version of a working pipeline, organised so that the logic, configuration and outputs are easy to reproduce and share.

## Workflow Summary

Inputs:

- A Cell Ranger multi CSV that points to your 10x FASTQ directories
- A gene-expression reference FASTA for souporcell
- A souporcell container image
- An arcasHLA environment with `arcasHLA` and `samtools`
- The `subset-bam` executable
- Optionally, a donor HLA annotation CSV to label clusters with real donor IDs

Outputs:

- Cell Ranger run directory under `results/cellranger/<run_id>/`
- Souporcell clusters under `results/souporcell/<run_id>/clusters.tsv`
- Donor barcode lists, donor BAMs, extracted FASTQs, and genotype JSONs under `results/arcasHLA/<run_id>/`
- HLA summary tables under `results/reports/<run_id>/`
- Donor-labeled demultiplexed outputs under `results/demultiplexed/<run_id>/`

If you provide a reference HLA annotation table, the workflow matches each souporcell cluster to the best donor profile using the configured HLA genes. If you leave that table blank, the workflow still runs and publishes outputs using the souporcell cluster names such as `donor0`, `donor1`, and so on.

## Repository Layout

```text
scHLA-demux/
├── config/
│   ├── config.yaml
│   ├── cellranger_multi.example.csv
│   └── reference_hla.example.csv
├── workflow/
│   ├── Snakefile
│   └── scripts/
│       ├── build_donor_assignments.py
│       ├── genotype_json_to_csv.py
│       └── publish_demultiplexed.py
└── README.md
```

## Requirements

The workflow expects these tools to be available on the machine where you run Snakemake:

- `snakemake`
- `cellranger`
- `apptainer` or `singularity` with a working souporcell container image
- `subset-bam`
- An arcasHLA environment whose `bin/` directory contains both `arcasHLA` and `samtools`

The helper scripts bundled in this repository only use the Python standard library.

## Input Files

### 1. Cell Ranger Multi CSV

Edit `config/cellranger_multi.example.csv` to point to your FASTQ directories, then update `inputs.cellranger_multi_csv` in `config/config.yaml`.

The workflow expects `cellranger multi`, the setup it was developed around. Only the Cell Ranger gene-expression BAM and barcode matrix are used downstream by souporcell and arcasHLA. If your experiment also includes VDJ, you can extend the CSV accordingly.

### 2. Optional Reference HLA Annotation CSV

If you want to map souporcell clusters back to known donor identities, create a CSV like `config/reference_hla.example.csv` and set `inputs.reference_hla_csv` in `config/config.yaml`.

Expected format:

- One row per donor
- A `donor` column
- One or more `*_alleles` columns such as `A_alleles`, `B_alleles`, `C_alleles`
- Multiple alleles for one gene separated by `|`

Example:

```csv
donor,A_alleles,B_alleles,C_alleles
Donor_01,A*02:01 | A*24:02,B*07:02 | B*44:03,C*05:01 | C*07:02
Donor_02,A*01:01 | A*03:01,B*08:01 | B*35:01,C*04:01 | C*07:01
```

## Configuration

Main settings live in `config/config.yaml`.

Important fields:

- `cellranger.run_id`: folder name created by Cell Ranger
- `cellranger.count_sample`: sample folder name under `outs/per_sample_outs/`
- `souporcell.donors`: expected donor count for souporcell
- `arcashla.env_bin`: `bin/` directory of the environment that contains `arcasHLA` and `samtools`
- `arcashla.subset_bam`: path to `subset-bam`
- `matching.genes`: HLA genes used to compare inferred cluster genotypes against your reference donor table

The default output structure is generated automatically under `results/`. You can keep the defaults, or override the souporcell and arcasHLA output directories in the YAML if you prefer a different layout.

## Run

From the repository root:

```bash
snakemake -s workflow/Snakefile --cores 42 --rerun-incomplete --printshellcmds
```

If you prefer to point to another config file:

```bash
snakemake -s workflow/Snakefile --configfile path/to/config.yaml --cores 42
```

## Final Outputs

The main user-facing outputs are:

- `results/reports/<run_id>/souporcell_cluster_hla_summary.csv`
  - One row per souporcell singlet cluster with the inferred HLA calls from arcasHLA
- `results/reports/<run_id>/cluster_to_donor_assignments.csv`
  - Cluster-to-donor assignment table
  - `match_status` tells you whether a cluster was confidently matched, ambiguous, unmatched, or published without a reference annotation
- `results/reports/<run_id>/donor_match_scores.csv`
  - Pairwise HLA match scores between each souporcell cluster and each reference donor
- `results/demultiplexed/<run_id>/<published_label>/`
  - `barcodes.txt`
  - `sample_alignments.bam`
  - `sample_alignments.bam.bai`
  - `reads.extracted.fq.gz`
  - `hla.genotype.json`
  - `hla.genes.json`
  - `hla.genotype.log`
  - `hla.alignment.p`

Those donor-specific folders are the “demultiplexed data ready to go” layer of the workflow. If a reference donor table is provided and the match is unique, the folder name is the donor ID from that table. Otherwise the workflow falls back to the souporcell cluster label.

## Notes

- The workflow reads only souporcell singlets when building donor-specific outputs.
- HLA matching defaults to `A,B,C`, because those genes are usually the most useful for donor identification, but you can change that in `matching.genes`.
- Matching is done at two-field HLA resolution by default so small formatting differences between reference annotations and arcasHLA output are less likely to block a valid match.
- If you already have Cell Ranger outputs in the expected location, the workflow reuses them instead of rerunning `cellranger multi`.

## Acknowledgements

scHLA-demux orchestrates existing tools; please cite them if you use it:
[Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger),
[souporcell](https://github.com/wheaton5/souporcell) (Heaton et al., *Nature Methods* 2020) and
[arcasHLA](https://github.com/RabadanLab/arcasHLA) (Orenbuch et al., *Bioinformatics* 2020).

## License

[MIT](LICENSE) © Roger Casals-Franch
