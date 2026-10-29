# scHLA-demux

`scHLA-demux` takes the FASTQ files of a 10x single-cell pool that mixes cells from several patients, splits the cells by patient, and tells you which patient each group of cells comes from, using their HLA type.

You give it:

- the FASTQ folders of each pool (gene expression, and optionally TCR),
- how many patients are in each pool,
- optionally, a table with each patient's known HLA alleles.

You get, for each pool:

- one folder per patient with that patient's cell barcodes, BAM, HLA calls and (if you sequenced TCR) TCR contigs,
- a table saying which patient every cell barcode belongs to, ready to add to Seurat or Scanpy,
- the HLA alleles found for each patient, in the same layout as the table you gave.

## How it works

```text
FASTQs (GEX [+ TCR])
   │
   ▼  cellranger multi         count matrix, BAM, TCR contigs
   ▼  souporcell               groups cells into N genotype clusters (N = patients in the pool)
   ▼  arcasHLA                 HLA alleles of each cluster
   ▼  HLA matching             cluster ↔ patient, using your HLA table
   ▼  publish                  one folder per patient + cell_to_donor.csv
```

Souporcell finds the patients from the SNPs in the RNA reads, but it only calls them `0, 1, 2, ...`. arcasHLA then genotypes the HLA genes of each cluster, and the matching step compares these alleles with your HLA table to put the real patient name on each cluster. Without an HLA table, everything still runs, and the patients keep the names `donor0`, `donor1`, ...

## Quick start

All commands are run from the `scHLA-demux` folder.

### 1. List your pools in `config/samples.csv`

One row per 10x pool (one Cell Ranger run):

```csv
pool,gex_fastqs,tcr_fastqs,n_donors,hla_csv
Pool1,/data/run1/L1216,/data/run1/L1213,4,config/hla_pool1.csv
Pool2,/data/run2/L1300,,3,
```

Here Pool1 has GEX and TCR, 4 patients and an HLA table. Pool2 is GEX only, with 3 patients and no HLA table. See [The sample sheet](#the-sample-sheet) for every column.

### 2. (Optional) Write an HLA table for each pool

```csv
donor,A_1,A_2,B_1,B_2,C_1,C_2
Patient_01,A*02:01,A*24:02,B*07:02,B*44:03,C*05:01,C*07:02
Patient_02,A*01:01,A*03:01,B*08:01,B*35:01,C*04:01,C*07:01
Patient_03,A*11:01,A*11:01,B*15:01,B*40:01,,
```

See [The donor HLA table](#the-donor-hla-table). `config/hla_example.csv` is a template.

### 3. Check, then run

```bash
# check every input and show what will run (small job, takes seconds)
sbatch --cpus-per-task=1 --mem=4G --partition=normal1 run_slurm.sh --dry-run

# run everything
sbatch run_slurm.sh
```

Logs are written to `logs/scHLA-demux_<jobid>.out` and `.err`. Before anything starts, the pipeline checks the sample sheet, the FASTQ folders, the HLA tables and the tool paths. Any problems are listed together in the `.err` log, so you can fix them all at once.

Results are written to `results/<pool>/`. Open `results/qc_report.html` first (see [Checking that everything ran well](#checking-that-everything-ran-well)). Then go to `results/<pool>/demultiplexed/`.

## The sample sheet

`config/samples.csv` has one row per pool. Comma-, semicolon- or tab-separated files all work, so a sheet saved from Excel is fine.

| Column | Required | What to put |
|---|---|---|
| `pool` | yes | A short name for the pool, used for the results folder. Letters, digits, `_`, `-`, `.` |
| `gex_fastqs` | yes | Folder with the gene-expression FASTQs (`NAME_S1_L001_R1_001.fastq.gz`, ...) |
| `tcr_fastqs` | no | Folder with the TCR (V(D)J-T) FASTQs. Leave empty for a GEX-only pool |
| `n_donors` | yes* | Number of patients mixed in the pool. *Can be left empty when `hla_csv` is given: the number of rows in the HLA table is used |
| `hla_csv` | no | The HLA table (`.csv`, `.tsv` or `.xlsx`) of the patients in this pool, or a table shared by several pools. Leave empty to skip matching |
| `gex_fastq_id`, `tcr_fastq_id` | no | Only needed when one folder holds FASTQs of several samples: the file name part before `_S1_...` |
| `expect_cells` | no | Expected number of cells for Cell Ranger (default: let Cell Ranger estimate) |
| `multi_csv` | no | Advanced: your own Cell Ranger multi CSV (e.g. for CITE-seq). The FASTQ columns are then ignored |

Notes:

- The pipeline writes the Cell Ranger multi CSV for you (`results/<pool>/cellranger_multi.csv`) and finds the FASTQ sample name from the file names.
- `n_donors` is the number souporcell looks for. Put the real number of patients in the pool, even if some are missing from your HLA table.
- Several pools in one sheet run in the same SLURM job, one after the other or side by side, as the cores allow.

## The donor HLA table

One row per patient, two columns per HLA gene (one allele in each):

```csv
donor,A_1,A_2,B_1,B_2,C_1,C_2
Patient_01,A*02:01,A*24:02,B*07:02,B*44:03,C*05:01,C*07:02
Patient_03,A*11:01,A*11:01,B*15:01,B*40:01,,
```

Rules:

- **`donor` column**: the patient name. This becomes the folder name and the label of the cells. `patient` also works as the column name.
- **Gene columns**: `<gene>_1` and `<gene>_2`, e.g. `A_1, A_2, B_1, B_2, C_1, C_2`. Matching uses A, B and C by default. You can add other genes such as `DRB1_1, DRB1_2`, and use them by adding them to `matching.genes` in `config/config.yaml`.
- **Homozygous**: write the same allele in both columns (`A*11:01, A*11:01`).
- **Unknown**: leave the cell empty (or write `NA`). That allele is just not compared.
- **Allele spelling**: `A*02:01`, `A*02:01:01`, `HLA-A*02:01` and `02:01` are all accepted. Matching is done at two-field resolution (`A*02:01`), so it doesn't matter how many fields you give.
- **Other columns** (e.g. `notes`) are allowed and ignored. The pipeline prints a warning so you can spot a typo such as `A1`.

Tip for Excel users: always write the gene prefix (`A*02:01`, not `02:01`). Otherwise Excel may turn the value into a time.

The pipeline writes the HLA alleles it finds in the same layout (`results/<pool>/reports/souporcell_cluster_hla_summary.csv`). You can open both tables side by side, or copy the output, rename the `donor` values and reuse it as an input.

**Other layouts that also work:**

- One column per gene, named after the gene, with both alleles in one cell. A single allele means homozygous:

  ```text
  patient      A                   B                   C         DRB1
  Patient_07   A*02:01, A*24:02    B*07:02, B*44:03    C*05:01   DRB1*04:01, DRB1*13:01
  ```

- `A_alleles` columns with alleles separated by `|`, the old layout.
- **File types**: `.csv`, `.tsv`, or an Excel `.xlsx` file. For Excel, the first sheet is read.

**One table for several pools**: you can point every pool to the same table. Each pool's clusters are matched against all patients in it, and the patients that belong to other pools are simply not used. In that case fill in `n_donors`, because otherwise it would be taken from the table's row count.

## Outputs

```text
results/
├── qc_report.html                     plots and checks for every pool: open this first
├── all_pools_donors.csv               one row per patient per pool: status, cell counts, HLA
├── all_pools_cell_to_donor.csv        every cell barcode of every pool, with its patient
└── <pool>/
    ├── demultiplexed/
    │   ├── cell_to_donor.csv          pool, barcode, donor, souporcell_status, souporcell_cluster
    │   ├── manifest.tsv               one row per patient, with paths to every file
    │   └── <patient>/
    │       ├── barcodes.txt                           this patient's cell barcodes
    │       ├── sample_alignments.bam(.bai)            this patient's reads
    │       ├── tcr_filtered_contig_annotations.csv    this patient's TCR contigs (TCR pools only)
    │       ├── hla.genotype.json                      arcasHLA call
    │       └── reads.extracted.fq.gz, hla.*           arcasHLA intermediate files
    ├── reports/
    │   ├── souporcell_cluster_hla_summary.csv   HLA alleles found per cluster (same layout as the input table)
    │   ├── cluster_to_donor_assignments.csv     which patient each cluster was matched to, and why
    │   └── donor_match_scores.csv               every cluster against every patient
    ├── cellranger_multi.csv                     the Cell Ranger config that was used
    ├── cellranger/                              full Cell Ranger output (web_summary.html, matrices, VDJ)
    ├── souporcell/                              full souporcell output
    └── arcasHLA/                                per-cluster intermediate files
```

In `cell_to_donor.csv`, the `donor` column holds the patient name for singlets, and `doublet` or `unassigned` for cells souporcell could not place.

### Checking that everything ran well

`results/qc_report.html` is a single file that opens in any browser, with nothing to install, and can be emailed. It starts with one line per pool. Each pool then gets:

- **Checks**, marked ✓ OK or ⚠ Check:
  - souporcell found the expected number of patients;
  - every patient has at least 200 cells;
  - doublets ≤ 30% and unassigned cells ≤ 10%;
  - at least 1,000 HLA reads per gene behind each HLA call;
  - every cluster matched a patient of the HLA table, and every patient of the table was found.
- **Cell Ranger numbers**: cells, genes and reads per cell, mapping, and TCR pairing.
- **Cells per patient**: singlets for each patient, plus doublets and unassigned cells.
- **HLA match grid**: for each cluster and each patient of the table, how many alleles they share (`5/5`). ✓ marks the name each cluster was given. A good run has one dark cell per row and per column.
- **TCR per patient**: the share of each patient's cells with a TCR contig (TCR pools only).
- **HLA calls**: the alleles found for each patient and the number of reads behind each gene.

⚠ Check means a result is unusual and worth a look, not that the run failed. Hover over a bar or grid cell to see its exact numbers. The thresholds are at the top of `workflow/scripts/qc_report.py`.

### `match_status`: how sure is the patient label?

| Status | Meaning | Folder name |
|---|---|---|
| `matched` | One patient clearly shares the most HLA alleles with this cluster | the patient name |
| `ambiguous` | Two or more patients tie for the best match (they share alleles) | `donorN` |
| `unmatched` | No allele shared with any patient in the table (e.g. patient missing from the table) | `donorN` |
| `conflicting_match` | Two clusters both matched best to the same patient | `donorN` |
| `no_reference_annotation` | No HLA table was given for this pool | `donorN` |

For anything other than `matched`, open `donor_match_scores.csv` to see how each cluster compares with each patient.

### Adding the patient labels to Seurat or Scanpy

The barcodes are the same as in Cell Ranger's filtered matrix.

```r
cells <- read.csv("results/Pool1/demultiplexed/cell_to_donor.csv")
obj$donor <- cells$donor[match(colnames(obj), cells$barcode)]
```

```python
import pandas as pd
cells = pd.read_csv("results/Pool1/demultiplexed/cell_to_donor.csv").set_index("barcode")
adata.obs["donor"] = cells["donor"].reindex(adata.obs_names).values
```

## Running, re-running and resuming

- **Resuming**: if the job stops (time limit, node failure), submit `sbatch run_slurm.sh` again. Finished steps are not redone.
- **Adding a pool**: add a row to `samples.csv` and submit again. Only the new pool runs.
- **Changing the HLA table**: edit it and submit again. Only the matching and publishing steps rerun.
- **Changing `n_donors`**: delete `results/<pool>/souporcell` and everything after it (`arcasHLA`, `reports`, `demultiplexed`), then submit again.
- **Re-running Cell Ranger**: an existing Cell Ranger output in `results/<pool>/cellranger` is always reused. Delete that folder to run it again (e.g. after changing the FASTQs).
- **Other config file**: `sbatch run_slurm.sh --configfile config/other.yaml`. Any extra arguments are passed to Snakemake.

## Configuration

Tool paths and settings live in `config/config.yaml`. You normally set them once per cluster:

| Setting | What it is |
|---|---|
| `samples` | Path of the sample sheet |
| `results_dir` | Where results are written |
| `cellranger.executable`, `gex_reference`, `vdj_reference` | Cell Ranger and its references (`vdj_reference` is only needed for TCR pools) |
| `cellranger.threads`, `localmem_gb` | Resources for Cell Ranger. Keep them within the SLURM request in `run_slurm.sh` |
| `souporcell.container` | Souporcell Apptainer/Singularity image |
| `souporcell.reference_fasta` | Genome FASTA. Empty = `<gex_reference>/fasta/genome.fa` |
| `arcashla.env_bin` | `bin/` folder of the conda env with `arcasHLA` and `samtools` |
| `arcashla.subset_bam` | 10x `subset-bam` executable |
| `arcashla.genes`, `population` | HLA genes to genotype and the arcasHLA population prior |
| `matching.genes` | HLA genes compared with your table (default `A,B,C`) |
| `matching.two_field_resolution` | Compare alleles at two-field resolution (default `true`) |

The SLURM request (partition, CPUs, memory, time) is set at the top of `run_slurm.sh`. The Snakemake used is set there too, or with `SNAKEMAKE=/path/to/snakemake sbatch run_slurm.sh`.

## Troubleshooting

- **"Please fix these problems before running scHLA-demux"**: the list under it says which line of which file is wrong.
- **"holds FASTQs of several samples"**: add a `gex_fastq_id` (or `tcr_fastq_id`) column with the sample name to use.
- **Cell Ranger says the pipestance is locked**: a previous run was killed. Delete `results/<pool>/cellranger/_lock` and submit again.
- **Most clusters are `unmatched`**: check that the HLA table belongs to this pool. Then compare `souporcell_cluster_hla_summary.csv` with your table. Low-coverage clusters can have incomplete HLA calls.
- **`ambiguous`**: the tied patients share alleles at the compared genes. Adding more genes (e.g. `DRB1`) to the table and to `matching.genes` can separate them.

## Requirements

- `snakemake` (7 or newer)
- `cellranger` (multi; tested with 9.0.1)
- `apptainer` or `singularity`, with the souporcell image
- `subset-bam` (10x Genomics)
- a conda env with `arcasHLA` and `samtools`

The helper scripts only use the Python standard library.

## Repository layout

```text
scHLA-demux/
├── run_slurm.sh                  SLURM launcher
├── LICENSE                       MIT
├── config/
│   ├── config.yaml               tool paths and settings
│   ├── samples.csv               the pools to run
│   ├── hla_example.csv           HLA table template
│   └── cellranger_multi.example.csv
└── workflow/
    ├── Snakefile
    └── scripts/
        ├── hla_table.py                reads/validates/writes HLA tables
        ├── make_multi_csv.py           writes the Cell Ranger multi CSV
        ├── genotype_json_to_csv.py     arcasHLA calls -> HLA table
        ├── build_donor_assignments.py  cluster <-> patient matching
        ├── publish_demultiplexed.py    per-patient folders, cell table, TCR split
        ├── summarize_pools.py          all-pools summary tables
        └── qc_report.py                QC report (HTML with inline SVG charts)
```

## Acknowledgements

scHLA-demux orchestrates existing tools; please cite them if you use it:
[Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger),
[souporcell](https://github.com/wheaton5/souporcell) (Heaton et al., *Nature Methods* 2020) and
[arcasHLA](https://github.com/RabadanLab/arcasHLA) (Orenbuch et al., *Bioinformatics* 2020).

## License

[MIT](LICENSE) © Roger Casals-Franch
