# SFARI project repository

This repository is a working analysis archive for the SFARI / transcription factor activation domain project. Most of the work lives in Jupyter notebooks, with a smaller set of reusable Python and R scripts for coordinate conversion, variant processing, plotting, and figure generation.

The goal of this README is to make handoff easier: if you need to find a notebook, script, figure, or output table, start here.

## Repository at a glance

### Top-level layout

| Path | What it contains |
| --- | --- |
| `/tmp/workspace/sanjanakotha/SFARI/notebooks` | Main analysis workspace. Most exploratory work, figure generation, conservation analyses, SFARI variant analyses, and activation-domain curation live here. |
| `/tmp/workspace/sanjanakotha/SFARI/output` | Exported tables, figures, alignments, predictor outputs, and per-project result folders. |
| `/tmp/workspace/sanjanakotha/SFARI/soto_analysis` | A more script-driven variant-analysis subproject, including helper scripts and notebooks for reproducing / extending the Soto-style gnomAD and patient-variant analyses. |
| `/tmp/workspace/sanjanakotha/SFARI/undergrad_projects` | Small side project / student-specific work. |
| `/tmp/workspace/sanjanakotha/SFARI/output_CC_15_ADs_upset_plot.png` | Standalone figure saved at the repo root. |

## Where to find things

### 1. Main notebook workspace: `notebooks/`

This is the core directory for the project. The notebook names are usually descriptive, so browsing by filename is often the fastest way to find prior work.

Common notebook themes in this folder:

- **SFARI-specific variant analyses**
  - `SFARI ADs.ipynb`
  - `SFARI ASD-Associated ADs.ipynb`
  - `SFARI TFs.ipynb`
  - `SFARI TF MTR.ipynb`
  - `SFARI variants and TF intersection .ipynb`
  - `SFARI variants intersection by TF.ipynb`
  - `Filtering SFARI TF variants.ipynb`
  - `SFARI GPF.ipynb`

- **Known AD list building / curation**
  - `Building List of All Known ADs.ipynb`
  - `Building List of All Known ADs- Redone Considering Isoforms.ipynb`
  - `Building GSL- Redone Considering Isoforms.ipynb`
  - `Building AD list - looking in canonical isoforms WITH ALERASOOL.ipynb`
  - `Adding alerasool ADs.ipynb`
  - `ADs per isoform.ipynb`
  - `Looking for non-canonical ADs.ipynb`

- **Predictor development / application to variants**
  - `AD Hunter on variants.ipynb`
  - `AD Hunter on Variants - SNP Changing V2.ipynb`
  - `ADPred results on variants.ipynb`
  - `PADDLE noSS results on variants.ipynb`
  - `Predictors - wt vs vars, continuous approach average tiles.ipynb`
  - `Predictors, continuous approach centered tiles.ipynb`
  - `AD tile classifier, activator TFs.ipynb`
  - `AD tile classifier, all TFs.ipynb`

- **Conservation / ortholog / phyloP analyses**
  - `100 Vertebrates.ipynb`
  - `100 Vertebrates All TFs.ipynb`
  - `100 Vertebrates, SFARI TFs.ipynb`
  - `Zoonomia.ipynb`
  - `Zoonomia All TFs.ipynb`
  - `Zoonomia SFARI TFs.ipynb`
  - `Zoonomia codon alignment.ipynb`
  - `Zoonomia alignment to residues.ipynb`
  - `Codon alignment, activator TFs.ipynb`
  - `Codon alignment, all TFs.ipynb`

- **Structure / disorder / constraint analyses**
  - `Disorder Predictions of Annotated ADs.ipynb`
  - `Disorder, Conservation, Constraint classifier.ipynb`
  - `Analyzing Activator TF DSSP.ipynb`
  - `Alphafold Activator TFs.ipynb`
  - `Gene Constraint Scores.ipynb`
  - `MTR All TFs.ipynb`

- **Project- or person-specific figure notebooks**
  - `Honors Thesis Figures.ipynb`
  - `Heatmaps for poster.ipynb`
  - `Cell conference figs.ipynb`
  - `Figures for Simons final report Dec 2025/`
  - `Summary figure per TF.ipynb`

- **Variant follow-up / experiment support**
  - `All AD annotations - Caitlin experiment.ipynb`
  - `Information on patient variants for Caitlin experiment - counts and alphamissense.ipynb`
  - `Other variants for Caitlin experiment - *.ipynb`
  - `Sending Caitlin Variants.ipynb`
  - `Splice AI Input.ipynb`
  - `ADs good for CRISPR deletion.ipynb`

### 2. Reusable helper code in `notebooks/`

A few Python files are reused across notebooks:

| File | Purpose |
| --- | --- |
| `/tmp/workspace/sanjanakotha/SFARI/notebooks/AD_predictor_tools.py` | Utilities for tiling protein sequences, computing amino-acid features, creating FASTA exports, and working with activation-domain predictor inputs. |
| `/tmp/workspace/sanjanakotha/SFARI/notebooks/AD_comparison_tools.py` | Functions for comparing predicted AD regions against gold-standard / reference sets and for overlap accounting. |
| `/tmp/workspace/sanjanakotha/SFARI/notebooks/PlottingTools.py` | Shared plotting helpers for histograms, heatmaps, and visual summaries used across notebook analyses. |
| `/tmp/workspace/sanjanakotha/SFARI/notebooks/summary_fig_plot_tools.py` | Higher-level figure helpers for summary visualizations, disorder/structure overlays, and publication-style plots. |
| `/tmp/workspace/sanjanakotha/SFARI/notebooks/uniprotBedTools.py` | Utilities for working with UniProt-linked coordinates and sequence regions inside notebooks. |
| `/tmp/workspace/sanjanakotha/SFARI/notebooks/ADpred_LambertTFs_helper.py` | Small helper module tied to ADpred/Lambert TF processing. |

If you are trying to rerun an older notebook and a function is missing, check these files first.

### 3. Script-driven analysis work: `soto_analysis/`

This subdirectory contains a more pipeline-like variant analysis workspace.

#### `soto_analysis/notebooks/`
These notebooks overlap conceptually with the main `notebooks/` folder, but are organized around the Soto / gnomAD / patient-variant analyses. Representative notebooks include:

- `Variant Enrichment Analysis.ipynb`
- `gnomAD analysis.ipynb`
- `gnomAD Variant Analysis.ipynb`
- `ClinVar Variant Analysis for 15.ipynb`
- `SFARI genome variants.ipynb`
- `SPARK iWES v2.ipynb`
- `SPARK iWES v3.ipynb`
- `Variants to fasta.ipynb`
- `Variants to predictor input.ipynb`
- `Rare Variants, AD vs TF.ipynb`
- `Rare Variants, AD vs DBD.ipynb`

#### `soto_analysis/scripts/`
This is where command-line helpers live.

| File | Purpose |
| --- | --- |
| `/tmp/workspace/sanjanakotha/SFARI/soto_analysis/scripts/classify_mutations.py` | Classifies SNVs as synonymous / nonsynonymous / nonsense by mapping genomic positions back to transcript codons. |
| `/tmp/workspace/sanjanakotha/SFARI/soto_analysis/scripts/coords_and_vars_to_seq.py` | Earlier sequence-mapping helper that converts genomic variants to codon/protein consequences. |
| `/tmp/workspace/sanjanakotha/SFARI/soto_analysis/scripts/get_protein_positions.py` | Maps nucleotide-level variant positions to protein positions. |
| `/tmp/workspace/sanjanakotha/SFARI/soto_analysis/scripts/get_full_mutations.py` | Builds mutation-level output files from variant inputs. |
| `/tmp/workspace/sanjanakotha/SFARI/soto_analysis/scripts/get_full_mutations_domains.py` | Variant extraction focused on domain-restricted outputs. |
| `/tmp/workspace/sanjanakotha/SFARI/soto_analysis/scripts/get_mutations_domains_snv_classified.py` | Domain-aware mutation processing plus SNV classification. |
| `/tmp/workspace/sanjanakotha/SFARI/soto_analysis/scripts/get_DBD_ED_coords.py` | Generates DNA-binding / effector-domain coordinate information. |
| `/tmp/workspace/sanjanakotha/SFARI/soto_analysis/scripts/gtf_to_bed.py` | Converts annotation features into BED-style coordinate files. |
| `/tmp/workspace/sanjanakotha/SFARI/soto_analysis/scripts/AD_coords_to_bed.R` | R helper for turning AD coordinates into BED format. |
| `/tmp/workspace/sanjanakotha/SFARI/soto_analysis/scripts/mapping_AA_coords_to_bed.R` | R helper for amino-acid-to-genome coordinate mapping. |

A lot of these scripts assume relative paths such as `../raw_files/` or `../outputs/`, so they are easiest to run from inside `soto_analysis/scripts/`.

### 4. Generated results: `output/`

This folder is the main landing spot for exported results from notebooks.

Important subdirectories / files:

| Path | What it appears to hold |
| --- | --- |
| `/tmp/workspace/sanjanakotha/SFARI/output/SFARI_known_ADs.csv` | Main SFARI known-AD table. |
| `/tmp/workspace/sanjanakotha/SFARI/output/SFARI_ADs_considering_isoforms_and_canonical.csv` | Curated SFARI AD list that includes isoform/canonical considerations. |
| `/tmp/workspace/sanjanakotha/SFARI/output/known_ADs_considering_isoforms_and_canonical*.csv` | Versions of the broader known-AD catalog used in later analyses. |
| `/tmp/workspace/sanjanakotha/SFARI/output/caitlin_experiment/` | Outputs prepared for the Caitlin-focused experiment: counts, ClinVar/gnomAD/COSMIC tables, AlphaMissense summaries, variant FASTAs, and figure exports. |
| `/tmp/workspace/sanjanakotha/SFARI/output/cc_TF_summary_vis/` | Per-transcription-factor summary figures (`.pdf` and `.png`). |
| `/tmp/workspace/sanjanakotha/SFARI/output/disorder_pred_figs/` | Disorder prediction figures for specific ADs / proteins. |
| `/tmp/workspace/sanjanakotha/SFARI/output/ensembl_orthologs/` | Ortholog FASTAs, MAFFT alignments, AD FASTAs, and percent-identity outputs. |
| `/tmp/workspace/sanjanakotha/SFARI/output/100_vertebrates_PhyloP_traces/` | 100-vertebrate conservation trace outputs. |
| `/tmp/workspace/sanjanakotha/SFARI/output/zoonomia_PhyloP_traces/` | Zoonomia conservation traces plus histogram/metapredict subfolders. |
| `/tmp/workspace/sanjanakotha/SFARI/output/variant_predictors/` | Predictor outputs on variants (for example `PADDLE_on_variants.csv`). |
| `/tmp/workspace/sanjanakotha/SFARI/output/nardini_output/` | Large collection of NARDINI PNG outputs. |
| `/tmp/workspace/sanjanakotha/SFARI/output/predictions/` | Saved prediction outputs. |

There are also many one-off CSVs, BED files, FASTA files, zipped outputs, and standalone figures directly under `output/`.

### 5. Small side area: `undergrad_projects/`

This directory currently contains a small amount of student-specific work (for example `test_CAMTA2_O94983_AD_285-468`). It does not look like the main project hub.

## Practical handoff notes

### How to navigate quickly

If you need to answer a question, these are good first stops:

- **Where was the SFARI TF / ASD analysis done?** Start in `notebooks/` with files beginning `SFARI ...`
- **Where was the known AD list assembled?** Start with the `Building ... ADs...` notebooks in `notebooks/`
- **Where were variant consequence / mutation mapping scripts kept?** Look in `soto_analysis/scripts/`
- **Where are exported figures and tables?** Look in `output/`
- **Where are per-gene summary figures?** Look in `output/cc_TF_summary_vis/`
- **Where are Caitlin-related exports?** Look in `output/caitlin_experiment/`
- **Where are conservation / ortholog outputs?** Look in `output/ensembl_orthologs/`, `output/100_vertebrates_PhyloP_traces/`, and `output/zoonomia_PhyloP_traces/`

### Data that is **not** in the repo

Some large or sensitive inputs are intentionally not tracked. In particular:

- `/tmp/workspace/sanjanakotha/SFARI/data/` is gitignored
- `/tmp/workspace/sanjanakotha/SFARI/soto_analysis/raw_files*` is gitignored
- `/tmp/workspace/sanjanakotha/SFARI/soto_analysis/output*` is gitignored

That means many notebooks and scripts may expect local files that are **not** present in this clone. If something fails because an input file is missing, check the notebook cells / script arguments for the expected filename and ask someone in the lab for the original data location.

### Environment / dependencies

There is currently **no tracked environment file** (`requirements.txt`, `environment.yml`, `pyproject.toml`, etc.) in the repository.

Based on the helper scripts, the main Python dependencies appear to include:

- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `scipy`
- `statsmodels`
- `biopython`
- `protfasta`
- `tqdm`
- `metapredict`

There is also some R-based work (`.R`, `.Rmd`) for coordinate mapping.

### Re-running old work

- Most analyses were run interactively in notebooks rather than as a single clean pipeline.
- File paths are often relative, so it helps to open notebooks from their original directory.
- Some notebooks are exploratory or iterative; if there are multiple similarly named notebooks, the newer / “redone” / “V2” / “considering isoforms” versions are usually the better starting points.
- `.ipynb_checkpoints/` directories are just notebook autosaves and can usually be ignored.

## Suggested starting points for a new lab member

1. Read this README once.
2. Skim `notebooks/` filenames to see the major analysis areas.
3. Open the specific notebook family relevant to your question:
   - `SFARI ...` for project-specific variant analyses
   - `Building ... ADs ...` for activation-domain catalog curation
   - `Zoonomia ...` or `100 Vertebrates ...` for conservation work
   - `AD Hunter ...`, `PADDLE ...`, or `ADPred ...` for predictor analyses
4. Use `output/` to find the final exported table or figure associated with that notebook.
5. If you need a scripted step rather than a notebook, check `soto_analysis/scripts/`.

## If I had to point you to just a few places

If you are totally new and only want the most useful folders:

- `/tmp/workspace/sanjanakotha/SFARI/notebooks`
- `/tmp/workspace/sanjanakotha/SFARI/output`
- `/tmp/workspace/sanjanakotha/SFARI/soto_analysis/scripts`

Those three locations cover most of the actual analysis history, exported results, and reusable processing code.
