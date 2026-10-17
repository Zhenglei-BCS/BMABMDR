# Simulation study: BMD10 / BMD20 / BMD50 – separate model averaging vs. one best-fitting model

This folder contains a self-contained simulation study (not part of the
installed package; it is excluded via `.Rbuildignore`) that compares two ways
of obtaining benchmark doses for several benchmark responses (BMR, `q` =
relative change in mean/median response) from **the same data** with BMABMDR:

* **Separate model averaging** – the way the package is designed to be used:
  `PREP_DATA_N/LN(q = 0.1)` + `full.laplace_MA()`, then again with
  `q = 0.2`, then with `q = 0.5` (three independent MA analyses).
* **One best-fitting model, all BMRs at once** – run the analysis once, take
  the model with the largest posterior weight and derive BMD10, BMD20 and
  BMD50 from that single fitted curve.

See [`PROGRESS.md`](PROGRESS.md) for a step-by-step work log (useful to resume
work), and the "Results" section below for the outcome of the run included in
this folder.

## 1. Feasibility

| Question | Finding |
|---|---|
| Can the package do "separate MA per BMR"? | Yes, natively: `q` is an argument of `PREP_DATA_N/LN`. One `full.laplace_MA` fit with all 16 continuous models takes ~1.5–4 s, so three fits per data set are cheap. Bridge sampling (`sampling_MA`) is far too slow for thousands of data sets and is not used. |
| Can the package do "one model, several BMRs"? | Not directly – every model is **parameterised with the BMD for the chosen `q` as a parameter** (`par2`; see `R/DRM.R`), and `q` also enters the prior constraint on the fold change (`L = 1+q+0.01` increasing, `U = 1-q-0.01` decreasing in `R/fun_Data.R`). However, `full.laplace_MA` returns the Laplace posterior draws of every model (`parsN`, `parsLN`; columns `p1..p4` = background, BMD, fold change, log d), and the exported `DRM.*` functions give the full fitted curve per draw. The BMD for any other BMR is then obtained by numerical inversion of each draw's curve (`bmd_from_draws()`), giving a full posterior – and thus BMDL/BMDU – for BMD20 and BMD50 from the fit at BMD10. |
| Is the inversion correct? | Checked: inverting each draw at the BMR used in the fit reproduces the fitted BMD draws to ~1e-11 for all 16 models, for increasing and decreasing data. |
| Are there BMRs that a single fit cannot provide? | Yes: a fit at `q = 0.1` only guarantees a fold change > 1.11 (increasing) / < 0.89 (decreasing). Posterior draws whose asymptote does not reach `q = 0.5` give BMD50 = ∞. This is reported as `p_inf` (share of non-estimable draws); infinite values are kept in the quantiles (so BMDU can be `Inf`). |
| Compute budget | ≈ 16 s per data set for all arms (3 MA fits + inversion) on one core. Full design (20 scenarios × 200 reps = 4000 data sets) ≈ 18 core-hours, i.e. ~4.5 h on 4 cores; the 30-replicate run included here takes ~50 min on 4 cores. |

Conclusion: the study is feasible with the Laplace approach and requires no
change to the package itself.

## 2. Design

### Estimators (arms)

For every simulated data set, with BMRs `q ∈ {0.1, 0.2, 0.5}`:

| Arm | Fits | Description |
|---|---|---|
| `MA_sep` | 3 | Separate `full.laplace_MA` per `q`; MA BMDL/BMD/BMDU (5%, 50%, 95%). |
| `Best_joint` | 1 | Fit at the reference BMR (`q_ref = 0.1`); model with the largest weight; BMD10/20/50 by inverting its posterior draws. |
| `MA_joint` | 1 (same) | Same single fit; BMDq derived for each model's draws and mixed according to the Laplace weights (model averaging of the derived BMDs). |
| `Best_sep` | 3 (same as `MA_sep`) | In each separate fit, the best model's own BMD posterior. |

`MA_sep` vs `Best_joint` is the comparison asked for. The two additional arms
make the comparison interpretable as a 2×2 factorial:
*averaging vs selection* (MA vs Best) × *one fit per BMR vs one fit for all
BMRs* (sep vs joint). At `q = q_ref` the joint arms coincide with the separate
ones (up to Monte-Carlo noise of the draws).

### Data-generating scenarios

Dose design `0, 6.25, 12.5, 25, 50, 100`, `n = 10` per group, background
`a = 10`. True curves (`true_curve()`):

| id | curve | parameters | true BMD10 / 20 / 50 |
|---|---|---|---|
| `E4_inc` | exponential | b=0.0144, c=3, d=1 | 3.56 / 7.32 / 19.98 |
| `H4_inc` | Hill | b=30, c=3, d=2 | 6.88 / 10.00 / 17.32 |
| `LN4_inc` | log-normal (probit in log dose) | b=40, c=3, d=1.5 | 13.36 / 17.02 / 25.51 |
| `LIN_inc` | linear, no plateau (misspecified for all models) | slope 0.02 | 5 / 10 / 25 |
| `H4_dec` | Hill, decreasing | b=30, c=0.3, d=2 | 12.25 / 18.97 / 47.43 |

crossed with error distribution (`normal`: constant SD = cv·a; `lognormal`:
SD of log = cv) and `cv ∈ {0.10, 0.20}` → 20 scenarios. True BMDs follow the
BMABMDR definition (mu(BMD) = mu(0)(1 ± q); for lognormal data the median),
computed by root finding (`true_bmd()`). Data are simulated per animal and
passed as dose-group summary statistics (mean, SD, n) to both
`PREP_DATA_N` and `PREP_DATA_LN`, with default priors and all 16 models.

### Performance measures (`summarise_performance()`)

Per scenario × BMR × arm: failure rate; median and mean of
log(BMD/true BMD) (bias) and RMSE on log scale; coverage of the 90 % credible
interval with Monte-Carlo SE; share with BMDL ≤ true BMD (protectiveness);
median log10(BMDU/BMDL) (precision); share of non-estimable draws; run time.
Per data set and arm: whether BMD10 ≤ BMD20 ≤ BMD50 and BMDL10 ≤ BMDL20 ≤
BMDL50 (internal consistency; `ordering_check()`). For `Best_sep`: whether
the selected model is the same at all three BMRs.

## 3. How to run

Requirements: BMABMDR installed (compiles the Stan models, ~20 min), plus
`ggplot2`. From the repository root:

```sh
# pilot: 2 scenarios x 3 replicates (~30 s on 4 cores)
Rscript simulation/run_simulation.R
Rscript simulation/analyze_results.R simulation/results/pilot

# full design (200 replicates per scenario)
SIM_MODE=full SIM_NREP=200 SIM_CORES=8 SIM_OUT=simulation/results/full_nrep200 \
  Rscript simulation/run_simulation.R
Rscript simulation/analyze_results.R simulation/results/full_nrep200
```

Other settings: `SIM_NDRAWS` (Laplace draws per model, default 10000),
`SIM_QREF` (reference BMR for the single-fit arms, one of 0.1/0.2/0.5).
Scenarios are defined at the top of `run_simulation.R`.
Every replicate is stored as `<out>/raw/<scenario>_repNNN.rds` (git-ignored);
**an interrupted run is resumed by re-running the same command**. Seeds are
deterministic per scenario and replicate.

Files:

* `R/sim_functions.R` – truth, data generation, BMD inversion, estimators,
  performance measures
* `run_simulation.R` – scenario grid, parallel loop, checkpointing, collection
* `analyze_results.R` – tables (`*.csv`) and figures (`fig_*.png`)
* `results/<run>/` – settings, true BMDs, collected results and summaries

## 4. Results

RESULTS_PLACEHOLDER
