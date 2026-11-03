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
   ▼  arcasHLA                 aligns each cluster's HLA reads to every known allele
   ▼  HLA from reads           HLA type of each cluster, and how well each patient's typing fits it
   ▼  HLA matching             cluster ↔ patient, using your HLA table
   ▼  publish                  one folder per patient + cell_to_donor.csv
```

**Use 10x 5' gene expression data** (e.g. 5' Immune Profiling, with or without TCR). HLA alleles differ mostly in exons 2 and 3, near the 5' end of the genes, which 5' libraries read: 87-92% of the HLA-A reads cover exons 1-3 in two public 10x 5' pools (Jurkat + MM.1S, and the 10x 5' hashing example). 3' libraries read the other end, where the alleles hardly differ: 0.5% of the HLA-A reads cover exons 1-3 in the 3' test data below. Souporcell still splits the patients, but most clusters cannot be named from the HLA table. The pipeline reads the chemistry from Cell Ranger and warns about 3' data.

Souporcell finds the patients from the SNPs in the RNA reads, but it only calls them `0, 1, 2, ...`. arcasHLA then aligns the HLA reads of each cluster to every known HLA allele, and the matching step checks, for each patient of your HLA table, how well their typing explains each cluster's reads, to put the real patient name on each cluster. Without an HLA table, everything still runs, and the patients keep the names `donor0`, `donor1`, ...

Why not use arcasHLA's final call directly? On 10x data it can go wrong: in its last step a null allele (an `N` allele, e.g. `A*24:608N`) can take all of a gene's reads, so the call comes out homozygous and the second allele is lost, although the reads show both clearly. scHLA-demux therefore goes back to the reads arcasHLA aligned, leaves null alleles out, and finds the allele pair that explains the most reads. The arcasHLA calls are still written, for reference (`souporcell_cluster_hla_summary.csv`).

## Installation

Three commands. `setup.sh` installs every tool and reference into `resources/` (about 50 GB, about an hour), and `config/config.yaml` already points there:

```bash
git clone https://github.com/rogercasalsfr/scHLA-demux.git
cd scHLA-demux
./setup.sh
```

You need a Linux machine (Cell Ranger runs only on Linux) with `curl`, `git` and `apptainer` (or `singularity`; ask your system administrator if it is missing). conda is used if you have it; if not, `setup.sh` first installs its own copy (Miniforge) into `resources/`, with no admin rights needed.

**Cell Ranger is the one part you fetch yourself**: it is free, but 10x gives the download link only after you accept their license. `setup.sh` asks for it at the start:

1. Open the [Cell Ranger download page](https://www.10xgenomics.com/support/software/cell-ranger/downloads) and fill in the form.
2. The page shows a `curl` command. Copy it, paste it when `setup.sh` asks (the whole command or just the link), and leave `setup.sh` running.

The link expires after about a day and is tied to your 10x account, so never save it in a file of the repo. If you already have `cellranger-<version>.tar.gz` (or `.tar.xz`), put it in `resources/` before running `./setup.sh`. To give the link without the question (e.g. in a script): `./setup.sh --cellranger-url "<link>"`, keeping the double quotes (the link contains `&`).

What `setup.sh` installs:

| Tool | Version | From |
|---|---|---|
| Cell Ranger | the one you download (tested: 10.1.0) | 10x Genomics (needs the link above) |
| GEX reference | `refdata-gex-GRCh38-2024-A` | 10x Genomics |
| V(D)J reference | `refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0` | 10x Genomics |
| souporcell | 2.5 | Docker image `cumulusprod/souporcell:2.5` (if that fails, the Sylabs library image `wheaton5/souporcell`) |
| `subset-bam` | 1.1.0 | 10x Genomics GitHub release |
| Snakemake | 9.27.1 | conda env `resources/envs/snakemake` (`envs/snakemake.yaml`) |
| arcasHLA, kallisto, samtools | 0.6.0, 0.44.0, 1.21 | conda env `resources/envs/arcashla` (`envs/arcashla.yaml`) |
| IMGT/HLA database | release 3.66.0 | the IMGT/HLA GitHub repository, fixed to that release |

Every version except Cell Ranger's is fixed in `setup.sh` and `envs/`, so a new installation gets the same tools and HLA database that the pipeline was tested with. Cell Ranger is the release you download from 10x; the pipeline reads the output layouts of Cell Ranger 7 to 10 and was tested with 10.1.0.

`setup.sh` can be run again at any time: finished parts are skipped, and an interrupted download or unpack is redone. `./setup.sh --check` only reports what is installed.

### Test your installation

One command runs the whole pipeline on a small public dataset and checks the result:

```bash
./test/run_test.sh      # about 8 GB download, then 2-4 hours; run it in tmux/screen or with nohup
```

The test data is the CellBench mixture of three lung cancer cell lines from three people (NCI-H1975, NCI-H2228, HCC827) in one 10x 3' v2 library, about 900 cells ([Tian et al. 2019](https://doi.org/10.1038/s41592-019-0425-8), GEO GSE111108). Their HLA types come from [Cellosaurus](https://www.cellosaurus.org), and the authors published which cell line each cell barcode comes from. Being 3' data, it tests the whole pipeline and the clustering, but cannot test the naming fully (see [How it works](#how-it-works)). The test passes when:

- souporcell finds the 3 cell lines: each cluster holds one line (≥ 95% of its cells), and at least 95% of the cells are in the cluster of their published line;
- the pipeline reports the 3' chemistry;
- every souporcell cluster gets a naming result (a name, `unmatched` or `ambiguous`);
- no cluster is named after the wrong cell line (with 3' data some stay unnamed, which is expected);
- with one cell line's HLA typing removed, still no cluster gets a wrong name.

It prints `TEST PASSED` or `TEST FAILED` with the checks, and the full results are in `test/results/` (start with `qc_report.html`). The test data has no TCR library, so the TCR steps are not covered.

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
# check every input and show what will run (takes seconds, no SLURM job)
./run_slurm.sh --dry-run

# run everything (one SLURM job on a whole node)
sbatch run_slurm.sh

# no SLURM: run it directly on the machine (uses all its CPUs; start it in tmux or screen)
./run_slurm.sh
```

Logs are written to `logs/scHLA-demux_<jobid>.out` and `.err`. Before anything starts, the pipeline checks the sample sheet, the FASTQ folders, the HLA tables and the tool paths. Any problems are listed together in the `.err` log, so you can fix them all at once.

Results are written to `results/<pool>/`. Open `results/qc_report.html` first (see [Checking that everything ran well](#checking-that-everything-ran-well)). Then go to `results/<pool>/demultiplexed/`.

## The sample sheet

`config/samples.csv` has one row per pool. Comma-, semicolon- or tab-separated files all work, so a sheet saved from Excel is fine.

| Column | Required | What to put |
|---|---|---|
| `pool` | yes | A short name for the pool, used for the results folder. Letters, digits, `_`, `-`, `.` |
| `gex_fastqs` | yes | The gene-expression FASTQs: their folder, or the folder plus the file name start when the folder holds several libraries (`/data/Pool1_GEX` for `/data/Pool1_GEX_R1.fastq.gz`). See [FASTQ names](#fastq-names) |
| `tcr_fastqs` | no | The TCR (V(D)J-T) FASTQs, in the same way. Leave empty for a GEX-only pool |
| `n_donors` | yes* | Number of patients mixed in the pool. *Can be left empty when `hla_csv` is given: the number of rows in the HLA table is used |
| `hla_csv` | no | The HLA table (`.csv`, `.tsv` or `.xlsx`) of the patients in this pool, or a table shared by several pools. Leave empty to skip matching |
| `gex_fastq_id`, `tcr_fastq_id` | no | Another way to choose the library in a folder with several: the file name part before `_S1_...` or `_R1` |
| `expect_cells` | no | Expected number of cells for Cell Ranger (default: let Cell Ranger estimate) |
| `manual_assignments` | no | A CSV naming clusters by hand (`cluster,donor`, e.g. `donor3,Patient_07`), for when HLA matching cannot decide. See [match_status](#match_status-how-sure-is-the-patient-label) |
| `multi_csv` | no | Advanced: your own Cell Ranger multi CSV (e.g. for CITE-seq). The FASTQ columns are then ignored |

Notes:

- The pipeline writes the Cell Ranger multi CSV for you (`results/<pool>/cellranger_multi.csv`) and finds the FASTQ sample name from the file names.
- `n_donors` is the number souporcell looks for. Put the real number of patients in the pool, even if some are missing from your HLA table.
- Several pools in one sheet run in the same SLURM job, one after the other or side by side, as the cores allow.

### FASTQ names

Both naming styles work:

- 10x / bcl2fastq names: `NAME_S1_L001_R1_001.fastq.gz`, `NAME_S1_L001_R2_001.fastq.gz` (several lanes are fine).
- Simple names, as many sequencing providers deliver them: `NAME_R1.fastq.gz`, `NAME_R2.fastq.gz` (also `NAME.R1.fq.gz`). The pipeline links them under 10x names in `results/<pool>/fastqs/` before Cell Ranger. The original files are not renamed, moved or copied.

R1 and R2 are required; index reads (I1, I2) are optional.

## The donor HLA table

One row per patient, two columns per HLA gene (one allele in each):

```csv
donor,A_1,A_2,B_1,B_2,C_1,C_2
Patient_01,A*02:01,A*24:02,B*07:02,B*44:03,C*05:01,C*07:02
Patient_03,A*11:01,A*11:01,B*15:01,B*40:01,,
```

Rules:

- **`donor` column**: the patient name. This becomes the folder name and the label of the cells. `patient` also works as the column name.
- **Gene columns**: `<gene>_1` and `<gene>_2` (or `<gene>1` and `<gene>2`), e.g. `A_1, A_2, B_1, B_2, C_1, C_2`. Matching uses A, B and C by default. You can add other genes such as `DRB1_1, DRB1_2`, and use them by adding them to `matching.genes` in `config/config.yaml`.
- **Homozygous**: write the same allele in both columns (`A*11:01, A*11:01`).
- **Unknown**: leave the cell empty (or write `NA`). That allele is just not compared.
- **Patients without HLA typing**: keep them in the table with empty HLA cells. If all the other patients of the pool are found, the cluster left over is named after them (`matched_by_elimination`).
- **Allele spelling**: `A*02:01`, `A*02:01:01`, `HLA-A*02:01` and `02:01` are all accepted.
- **P and G groups** (`B*39:01P`, `A*02:01:01G`, as many HLA labs report): the HLA reference has no group names, so the P or G is dropped and the group is used as the allele that names it (`B*39:01`, `A*02:01:01`). The pipeline lists these conversions at the start of the run.
- **Low-resolution typing works**: `A*03`, `A3`, or just `3` in an A column (also `Cw7` for C). Each allele is used at the resolution it has: `A*03` in the table allows any `A*03` allele, `A*03:01` only `A*03:01` alleles (the table is used at two fields at most). A gene with only one allele typed leaves the other one free.
- **`pool` column** (optional): which pool each patient is in. With it, one table can hold the patients of every pool: each pool is only matched against its own patients, and `n_donors` can be left empty.
- **Other columns** (age, diagnosis, notes...) are allowed and ignored, so a patient sheet with clinical data works as is. The pipeline only warns about a column that looks like a mistyped HLA column, such as `A3`.

Tip for Excel users: always write the gene prefix (`A*02:01`, not `02:01`). Otherwise Excel may turn the value into a time.

The pipeline writes the HLA types it finds in the same layout (`results/<pool>/reports/cluster_hla_from_reads.csv`, at one field). You can open both tables side by side, or copy the output, rename the `donor` values and reuse it as an input.

**Other layouts that also work:**

- One column per gene, named after the gene, with both alleles in one cell. A single allele means homozygous:

  ```text
  patient      A                   B                   C         DRB1
  Patient_07   A*02:01, A*24:02    B*07:02, B*44:03    C*05:01   DRB1*04:01, DRB1*13:01
  ```

- `A_alleles` columns with alleles separated by `|`, the old layout.
- **File types**: `.csv`, `.tsv`, or an Excel `.xlsx` file. For Excel, the first sheet is read.

**One table for several pools**: add a `pool` column (see above), or point every pool to the same table without it. Without a `pool` column, each pool's clusters are matched against all patients in the table, and the patients that belong to other pools are simply not used. In that case fill in `n_donors`, because otherwise it would be taken from the table's row count.

## Outputs

```text
results/
├── qc_report.html                     plots and checks for every pool: open this first
├── all_pools_donors.csv               one row per patient per pool: status, cell counts, HLA
├── all_pools_cell_to_donor.csv        every cell barcode of every pool, with its patient
├── run_info/<date_time>/              copy of the config, sample sheet and HLA tables, and the tool versions, for each run
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
    │   ├── cluster_hla_from_reads.csv           HLA type of each cluster from its reads (same layout as the input table)
    │   ├── hla_read_evidence.csv                per cluster and gene: reads, call, support of each allele, closest alternative
    │   ├── donor_read_fit.csv                   per cluster, patient and gene: how well the patient's typing fits the reads
    │   ├── souporcell_cluster_hla_summary.csv   arcasHLA's own calls, for reference
    │   ├── cluster_to_donor_assignments.csv     which patient each cluster was matched to, and why
    │   └── donor_match_scores.csv               every cluster against every patient
    ├── logs/                                    one log per step (kept when a step fails)
    ├── fastqs/                                  10x-named links to your FASTQs (only for simple names)
    ├── cellranger_multi.csv                     the Cell Ranger config that was used
    ├── cellranger/                              full Cell Ranger output (web_summary.html, matrices, VDJ)
    ├── souporcell/                              full souporcell output
    └── arcasHLA/                                per-cluster intermediate files
```

In `cell_to_donor.csv`, the `donor` column holds the patient name for singlets, and `doublet` or `unassigned` for cells souporcell could not place.

### Checking that everything ran well

`results/qc_report.html` is a single file that opens in any browser, with nothing to install, and can be emailed. It starts with one line per pool. Each pool then gets:

- **Checks**, marked ✓ OK or ⚠ Check:
  - the library is 5' gene expression (3' data: HLA naming mostly fails, see [How it works](#how-it-works));
  - souporcell found the expected number of patients;
  - every patient has at least 200 cells;
  - doublets ≤ 30% and unassigned cells ≤ 10%;
  - at least 1,000 HLA reads per gene behind each HLA call;
  - every cluster matched a patient of the HLA table, and every patient of the table was found;
  - the named patients' HLA typing fits their cluster's reads (fit ≥ 0.95), and every other naming fits clearly worse (margin ≥ 0.05).
- **Cell Ranger numbers**: cells, genes and reads per cell, mapping, and TCR pairing.
- **Cells per patient**: singlets for each patient, plus doublets and unassigned cells.
- **HLA match grid**: for each cluster and each patient of the table, how well the patient's typing fits the cluster's HLA reads (`1.00` = as well as the best possible HLA type; a wrong patient is usually below `0.9`). ✓ marks the name each cluster was given. A good run has one dark `1.00` cell per row and per column.
- **TCR per patient**: the share of each patient's cells with a TCR contig (TCR pools only).
- **HLA calls**: the HLA type found for each patient, the support of each allele (see below), the fit and margin of the naming, and the number of reads behind each gene.

⚠ Check means a result is unusual and worth a look, not that the run failed. Hover over a bar or grid cell to see its exact numbers. The thresholds are at the top of `workflow/scripts/qc_report.py`.

### `match_status`: how sure is the patient label?

Each cluster is scored against each patient by the **read fit**: the share of the cluster's HLA reads explained by the best allele pair the patient's typing allows, divided by the share explained by the best allele pair overall, averaged over the matching genes. `1.000` means the typing explains the reads as well as any HLA type can; a patient whose typing misses one of the cluster's alleles drops clearly (in the public Jurkat + MM.1S mixture, each cell line's typing fits its own cluster at 0.995-1.000 and the other cluster at 0.32-0.33). A typed allele with no reads costs nothing, so a person who carries the cluster's alleles plus another one can fit as well as the true patient: the naming then comes out `ambiguous`, not wrong.

Clusters are then named one-to-one: each patient names at most one cluster, and the pipeline picks the naming with the best total fit over the whole pool. So two clusters can never get the same patient, and a cluster that only resembles a patient (shares a common allele such as `A*02`) does not take the name of a patient that another cluster fits better. A cluster is only named after a patient whose fit is at least `matching.min_read_fit` (0.90).

| Status | Meaning | Folder name |
|---|---|---|
| `matched` | Named from the HLA comparison, and no other naming fits nearly as well | the patient name |
| `matched_by_elimination` | Nothing to compare (the patient has no HLA typing, or the cluster has too few HLA reads), but this patient is the only one left once the other clusters are named | the patient name |
| `manual` | Named in your `manual_assignments` file | the patient name |
| `ambiguous` | Two or more namings fit (almost) equally well: margin below `matching.min_margin`, e.g. patients who share their alleles at the typed resolution, two patients without typing, or one patient without typing and two clusters left without a match | `donorN` |
| `unmatched` | No patient of the table fits well enough (fit below `matching.min_read_fit`), or there are more clusters than patients | `donorN` |
| `no_reference_annotation` | No HLA table was given for this pool | `donorN` |

How clear each naming is, in `cluster_to_donor_assignments.csv` (and `all_pools_donors.csv`):

- `read_fit`: fit of the named patient (`read_fit_by_gene` per gene);
- `naming_margin`: how much lower the total fit is if this cluster is given to any other typed patient instead. In a clean run it is well above 0.1;
- `runner_up`: the next best patient for this cluster, with their fit.

Open `donor_match_scores.csv` to see how each cluster compares with each patient, and `hla_read_evidence.csv` for the HLA evidence of each cluster: the HLA type from the reads, the `allele_support_pct` of each allele (the % of the gene's reads lost when that allele is swapped for the best other one: small means a similar allele, such as `A*23` for `A*24`, fits almost as well), the closest alternative type, and the reads that fit only one allele group. A second allele that adds fewer than 2% of the reads is not called (the gene is reported homozygous): such reads are usually ambient RNA from the other patients of the pool.

The older matching, which compares arcasHLA's final calls allele by allele with the table, is still available with `matching.method: alleles`.

**Naming clusters by hand**: when HLA matching cannot decide (`ambiguous`), or you know a cluster from other data (sex genes, a known mutation...), write a small CSV and add it as `manual_assignments` in the sample sheet:

```csv
cluster,donor
donor3,Patient_07
```

The other clusters are then named around these, and only the matching and publishing steps rerun.

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

### Running one stage at a time

Every step keeps its output (nothing is deleted at the end), so you can stop after any stage, look at the results, and continue later from there:

```bash
sbatch run_slurm.sh --until cellranger_multi    # stop after Cell Ranger: check web_summary.html
sbatch run_slurm.sh --until run_souporcell      # stop after souporcell: check clusters.tsv
sbatch run_slurm.sh                             # continue to the end
```

To redo a single step, delete its output (or use `--forcerun <step>`) and submit again. The steps are, in order: `link_fastqs`, `make_multi_csv`, `cellranger_multi`, `run_souporcell`, `make_donor_barcodes`, `filter_donor_bam`, `extract_donor_fastq`, `genotype_donor`, `summarize_genotypes`, `hla_read_evidence`, `build_donor_assignments`, `publish_demultiplexed`, `summarize_pools`, `qc_report`.

### When a step fails

- Every step writes its own log to `results/<pool>/logs/<step>.log` (per cluster for the HLA steps, e.g. `donor2_arcashla_genotype.log`). The logs are kept when a step fails, and the SLURM `.err` log names the step and its log file.
- The other pool keeps running (`--keep-going`), so one broken pool does not stop the rest.
- Fix the problem and submit again: finished steps are kept, and the failed one restarts. Cell Ranger also continues its own run where it stopped.

## Configuration

Tool paths and settings live in `config/config.yaml`. The tool paths point to what `setup.sh` installed in `resources/`; change them only to use tools installed elsewhere:

| Setting | What it is |
|---|---|
| `samples` | Path of the sample sheet |
| `results_dir` | Where results are written |
| `cellranger.executable`, `gex_reference`, `vdj_reference` | Cell Ranger and its references (`vdj_reference` is only needed for TCR pools) |
| `cellranger.threads`, `localmem_gb` | Resources for Cell Ranger. `auto` (default) uses every CPU of the job and 90% of its memory |
| `souporcell.container` | Souporcell Apptainer/Singularity image |
| `souporcell.reference_fasta` | Genome FASTA. Empty = `<gex_reference>/fasta/genome.fa` |
| `arcashla.env_bin` | `bin/` folder of the conda env with `arcasHLA` and `samtools` |
| `arcashla.subset_bam` | 10x `subset-bam` executable |
| `arcashla.genes`, `population` | HLA genes to genotype and the arcasHLA population prior |
| `matching.method` | `reads` (default): score each patient's typing against the cluster's HLA reads. `alleles`: compare arcasHLA's calls with the table |
| `matching.genes` | HLA genes compared with your table (default `A,B,C`) |
| `matching.two_field_resolution` | Use the table at two-field resolution at most (default `true`) |
| `matching.min_read_fit`, `min_margin` | Lowest read fit for naming a cluster (default `0.90`), and the margin below which a naming is `ambiguous` (default `0.03`) |

`run_slurm.sh` asks SLURM for a whole node (all CPUs and memory) for up to 7 days. To change that, add `sbatch` options, which override the ones in the script: `sbatch --partition=highmem --time=2-00:00:00 run_slurm.sh`. The Snakemake used is set there too, or with `SNAKEMAKE=/path/to/snakemake sbatch run_slurm.sh`.

## Troubleshooting

- **"Please fix these problems before running scHLA-demux"**: the list under it says which line of which file is wrong.
- **"holds FASTQs of several samples"**: put the folder plus the file name start in `gex_fastqs` (e.g. `/data/Pool1_GEX` for `/data/Pool1_GEX_R1.fastq.gz`), or add a `gex_fastq_id` (or `tcr_fastq_id`) column with that name.
- **A run was killed** (out of memory, node reboot, time limit): submit again. Cell Ranger's leftover lock is removed automatically and it resumes where it stopped.
- **Most clusters are `unmatched`**: first check the chemistry line of the QC report: with 3' data this is expected (name the clusters with `manual_assignments` if you know them from other data). Otherwise check that the HLA table belongs to this pool. Then compare `cluster_hla_from_reads.csv` with your table, and look at `read_fit_by_gene` in `donor_match_scores.csv` to see which gene disagrees: a typing error in the table shows up as one gene with a low fit. Low-coverage clusters (few HLA reads, e.g. 3' data) give less certain calls.
- **`ambiguous`**: the tied patients share alleles at the compared genes. Adding more genes (e.g. `DRB1`) to the table and to `matching.genes` can separate them.

## Requirements

`setup.sh` installs everything (see [Installation](#installation)). The machine needs:

- `curl`, `git` and `apptainer` (or `singularity`); conda is optional (`setup.sh` installs Miniforge without it)
- SLURM for `sbatch run_slurm.sh`, or any Linux machine to run `./run_slurm.sh` directly
- about 50 GB for `resources/`, plus 50-100 GB per pool for the results

The helper scripts only use the Python standard library, except `hla_read_evidence.py` (numpy), which runs in the arcasHLA env that `setup.sh` installs.

## Repository layout

```text
scHLA-demux/
├── setup.sh                      installs the tools and references into resources/
├── run_slurm.sh                  SLURM launcher
├── LICENSE                       MIT
├── envs/                         conda envs installed by setup.sh
├── test/                         ./test/run_test.sh: full run on public data, checked against the ground truth
├── resources/                    tools and references (created by setup.sh, not in git)
├── config/
│   ├── config.yaml               tool paths and settings
│   ├── samples.csv               the pools to run
│   ├── hla_example.csv           HLA table template
│   └── cellranger_multi.example.csv
└── workflow/
    ├── Snakefile
    └── scripts/
        ├── hla_table.py                reads/validates/writes HLA tables
        ├── chemistry.py                10x chemistry from Cell Ranger (warns about 3' data)
        ├── make_multi_csv.py           writes the Cell Ranger multi CSV
        ├── genotype_json_to_csv.py     arcasHLA calls -> HLA table
        ├── hla_read_evidence.py        HLA type from the reads, and fit of each patient's typing
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
