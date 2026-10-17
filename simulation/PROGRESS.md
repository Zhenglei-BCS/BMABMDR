# Simulation study – progress log

Purpose: make it possible to resume this work at any time. Every step is
logged here (newest at the bottom). If work stopped unexpectedly, read the
"Next steps" section at the end first.

Study question: does estimating BMD10, BMD20 and BMD50 with **three separate
BMABMDR model-averaging fits** (one per BMR `q`) perform differently from
estimating all three **from one fit of the best model** (and, as extra
reference arms, from one MA fit, or from the best model of each separate fit)?

## Log

### Step 1 – Repository exploration (done)
- Main estimator: `full.laplace_MA(data.N, data.LN, ...)` (R/fun_Laplace.R).
  Bridge sampling (`sampling_MA`) is far too slow for a large simulation; the
  study uses the Laplace approach.
- The BMR `q` is fixed at data preparation (`PREP_DATA_N/LN(q = ...)`,
  R/fun_Data.R). The BMD is a *model parameter* (`par2`, on dose scale
  rescaled by max dose) and `q` enters:
  - the curve parameterisation (`DRM.*` in R/DRM.R: `b` is a function of BMD, q),
  - the prior/constraint on fold change: increasing `L = 1+q+0.01`,
    decreasing `U = 1-q-0.01`, i.e. the fit at q forces the asymptote to be
    beyond a q-relative change.
- BMD definition (continuous): Normal models: mu(BMD) = mu(0)(1 +/- q);
  Lognormal models: median(BMD) = median(0)(1 +/- q)
  (log-scale `DRM.*_LN*` + log(1 +/- q)).
- Output of `full.laplace_MA`: `MA` (BMDL/BMD/BMDU), `weights` (16 Laplace
  weights in `get_models('continuous')` order), `BMDMixture`, `parsN`/`parsLN`
  (per-model data.frames of Laplace posterior draws, columns p1..p4 = mu(0),
  BMD (rescaled), fold change, log(d)), `E4_N` etc. (per-model summaries; row
  `BMD` holds BMDL/BMD/BMDU on original scale), `max.dose`, `q`, `shift`.
- `DRM.<model><I|D>(par = p1..p4, x, q, shift)` returns the fitted mean
  (Normal) or log-median (Lognormal) curve, so a BMD for *any* other BMR can be
  obtained per posterior draw by numerical root-finding on the fitted curve.
- R was not installed in the sandbox; installed via apt (r-cran-*) + CRAN
  (`mc2d`, `gamlss`), then the package has to be compiled (47 Stan models).

### Step 2 – Environment (done)
- `apt-get install r-base-core r-base-dev r-cran-rstan r-cran-rstantools
  r-cran-bridgesampling r-cran-bayesplot r-cran-ggplot2 r-cran-ggpubr
  r-cran-dplyr r-cran-plyr r-cran-shiny r-cran-rcolorbrewer r-cran-mvtnorm
  r-cran-posterior r-cran-pracma r-cran-truncnorm r-cran-rcppparallel
  r-cran-rcppeigen r-cran-bh r-cran-stanheaders r-cran-testthat`
- CRAN is blocked in the sandbox: `gamlss.data`, `gamlss.dist`, `gamlss`,
  `mc2d` installed from the GitHub CRAN mirror
  (`https://codeload.github.com/cran/<pkg>/tar.gz/refs/heads/master`, then
  `R CMD INSTALL`).
- Package installed with `MAKEFLAGS=-j4 R CMD INSTALL .` (compiles Stan models).

### Step 3 – Design decisions (done)
Four estimators per simulated data set (BMRs q = 0.1, 0.2, 0.5):

| Arm          | # fits | How BMDq is obtained |
|--------------|--------|----------------------|
| `MA_sep`     | 3      | `full.laplace_MA` run separately with q = 0.1, 0.2, 0.5 (MA BMDL/BMD/BMDU) |
| `Best_sep`   | (same 3) | in each separate fit, the model with the largest weight, its own BMD posterior |
| `Best_joint` | 1      | fit at q_ref (default 0.1), best model; BMD10/20/50 from inverting its posterior curve draws |
| `MA_joint`   | (same 1) | fit at q_ref; BMDq derived for every model's draws, mixed by Laplace weights |

`Best_joint` vs `MA_sep` is the comparison asked for; the other two arms
disentangle *model averaging vs model selection* from *separate vs single fit*.
Code: `simulation/R/sim_functions.R` (truth, data generation, inversion
`bmd_from_draws()`, estimators, performance measures).

### Step 4 – Scripts written (done)
- `simulation/run_simulation.R`: scenario grid, parallel (`mclapply`) loop,
  **one .rds per scenario x replicate** in `<out>/raw/` (re-running resumes),
  collection into `<out>/all_results.rds`.
- `simulation/analyze_results.R`: performance tables + figures.
- Note: package install must be done from a copy because `configure` in the
  repo is not executable:
  `git archive HEAD | tar -x -C /tmp/bmabuild && chmod +x /tmp/bmabuild/configure*
   && MAKEFLAGS=-j4 R CMD INSTALL /tmp/bmabuild` (~15-20 min, 46 Stan models).
- True BMDs checked against closed forms (e.g. Hill d=2, b=30, c=3:
  BMD20 = 30*sqrt(0.1/0.9) = 10, BMD50 = 30/sqrt(3) = 17.32) – OK.

### Step 5 – Self-checks and pilot (done)
- `bmd_from_draws()` evaluated at q_new = q_ref reproduces the fitted BMD
  draws (`par2`) for all 16 models, increasing and decreasing data:
  max abs difference ~2.5e-11.
- One `full.laplace_MA` fit (16 models, ndraws = 5000-10000) takes ~1.5-4 s.
- Pilot (`Rscript simulation/run_simulation.R`, 2 scenarios x 3 reps):
  6 tasks in 0.5 min on 4 cores, all four arms produce estimates. Pilot
  output was not kept (it is regenerated in 30 s).
- `simulation/` added to `.Rbuildignore`; raw replicate files are git-ignored
  (`simulation/.gitignore`).

### Step 6 – Full-grid run with 30 replicates (started 2026-10-17 ~18:55 UTC)
Command:
`SIM_MODE=full SIM_NREP=30 SIM_CORES=4 SIM_OUT=simulation/results/full_nrep30 Rscript simulation/run_simulation.R`
then `Rscript simulation/analyze_results.R simulation/results/full_nrep30`.
If interrupted: simply re-run the same command; finished replicates in
`raw/` are skipped.
