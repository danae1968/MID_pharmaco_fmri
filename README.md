# Dopaminergic drug effects on reward processing: MID fMRI analysis

Analysis code for a within-subject, double-blind, placebo-controlled
pharmacological fMRI study. About 100 participants performed the
**monetary incentive delay (MID) task** in three sessions, after
**methylphenidate** (MPH, a dopamine/noradrenaline reuptake inhibitor),
**sulpiride** (SUL, a D2 receptor antagonist) or **placebo** (PBO). Baseline
striatal dopamine synthesis capacity (*ki*) was measured with [18F]-FDOPA
PET.

The code runs from per-subject GLMs through permutation-based group
inference (TFCE) to publication tables and figures, along with the
behavioural analysis and brain–behaviour plots.

**Tools:** MATLAB · SPM12 · TFCE toolbox · R (afex, emmeans, tidyverse, ggplot2) · Slurm

## Pipeline

| Step | Folder | Entry point | What it does |
|---|---|---|---|
| 1 | `01_first_level/` | `run_1stLevel_MID.m` | Task onsets from log files (`onsetsMID.m`); per-session GLM (canonical HRF, fMRIPrep confounds); task contrasts; per-subject drug-difference and session-average images |
| 2 | `02_behaviour/` | `MID_analysis.Rmd` | Outlier exclusion; Reward × Drug repeated-measures ANOVAs; *ki* moderation (ANCOVA); ex-Gaussian RT modelling; per-subject behavioural drug effects for step 3 |
| 3 | `03_group_level/` | `run_groupLevel_vici_MID.m` | One-sample t-tests on session-average, drug-comparison and single-drug images, optionally with a *ki* or behavioural covariate; significance masks |
| 4 | `04_tfce/` | `wrapper_tfce.m`, then `write_binary_tfce.m` | Threshold-free cluster enhancement (5000 permutations, FWE p < .05), restricted to task-responsive voxels for drug and covariate effects; one Slurm job per model |
| 5 | `05_cluster_tables/` | `wrapper_extract_clusters_MID.m`, then `combine_cluster_tables.m` | Cluster/peak tables with atlas labels and per-region breakdowns |
| 6 | `06_figures/` | `wrapper_slice_display.m` | Dual-coded slice figures (effect size as colour, \|t\| as opacity, significance as contour) |
| 7 | `07_brain_behaviour/` | `extract_neural_beh_link.m`, then `neural_beh_scatter.R` | Per-subject mean contrast in the principal parietal cluster of each behavioural covariate model, plotted against the behavioural drug effect (Figure 8B, S1B) |

Run the steps in order. Two points of order:

- In step 3, run the session-average models without a covariate first.
  Their p < .001 masks serve as the functional mask in steps 3–4.
- Knit `MID_analysis.Rmd` before step 3, which needs
  `MID_beh_covariates.csv`. Knit it again after step 5 to build the fMRI
  cluster tables.

Each entry-point script starts with a configuration block (analyses,
contrasts, covariates, thresholds).

## Requirements

- MATLAB (R2019a or later) with [SPM12](https://www.fil.ion.ucl.ac.uk/spm/)
  and the [TFCE toolbox](https://github.com/ChristianGaser/tfce)
- [slice_display](https://github.com/bramzandbelt/slice_display) and
  [panel](https://www.mathworks.com/matlabcentral/fileexchange/20003-panel)
  for step 6
- R ≥ 4.1 with `dplyr`, `tidyr`, `readr`, `ggplot2`, `afex`, `emmeans`,
  `rstatix`, `ggpubr`, `ggrain`, `kableExtra`, `effectsize`, `sjPlot`

In MATLAB, run `setup_paths.m` once per session to add all step folders to
the path.

## Data

The data are not included. Paths point to the project storage of the
Donders Institute (BIDS layout, fMRIPrep derivatives). The R scripts can
read a local copy with the same layout via the `MID_LOCAL_DATA`
environment variable.

## Author and credits

Danae Papadopetraki, README by Claude

The first- and group-level scripts started from a lab template by Ruben
van den Bosch, and `estimate_tfce_mid.m` is adapted from his TFCE
function. His own analysis code for the same project is available at
<https://doi.org/10.34973/bc23-mz79> (van den Bosch et al., 2022,
*Nature Communications*, <https://doi.org/10.1038/s41467-022-32679-1>).
