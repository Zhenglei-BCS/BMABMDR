# Bug: `full.laplaceQ_MA()` crashes with "missing value where TRUE/FALSE needed" when all models fail

## Summary

When fitting clustered (Beta-Binomial) quantal data with `full.laplaceQ_MA()`,
if **all 8 models** fail the Laplace optimisation (i.e. every `rstan::optimizing()`
call hits a non-positive-definite Hessian), the model-averaged posterior weights
are all `NA`. The downstream code then crashes at:

```r
if (maci[2]/maci[1] > 20) ...
#> Error: missing value where TRUE/FALSE needed
```

This is a separate issue from the `fpfit2$sigma.coefficients` bug (see
`bug_report_gamlss_try_error.md`), though the two often co-occur on sparse
clustered data.

## Reproducible example

```r
library(BMABMDR)

# Sparse simulated data: only 3 litters per dose
set.seed(42)
indiv <- data.frame(
  dose   = rep(c(0, 5, 15, 50, 100), each = 60),
  litter = rep(rep(1:3, each = 20), 5),
  y      = c(rbinom(60, 1, 0.05),
             rbinom(60, 1, 0.10),
             rbinom(60, 1, 0.20),
             rbinom(60, 1, 0.35),
             rbinom(60, 1, 0.60))
)

library(dplyr)
litter_dat <- indiv |>
  group_by(dose, litter) |>
  summarise(y = sum(y), n = n(), .groups = "drop") |>
  select(dose, y, n) |>
  as.data.frame()

# 15 rows — too sparse for Beta-Binomial Laplace
data_QL <- PREP_DATA_QA(litter_dat, sumstats = TRUE, q = 0.1, cluster = TRUE)
fit <- full.laplaceQ_MA(data.Q = data_QL, prior.weights = rep(1, 8))
#> [many "Error in chol.default(-H)" messages]
#> Error in `if (maci[2]/maci[1] > 20) ...`:
#> ! missing value where TRUE/FALSE needed
```

## Root cause

In `full.laplaceQ_MA()` (file `R/fun_Laplace.R`), after attempting all 8 models,
the code computes model-averaged CrI bounds (`maci`). When every model's Laplace
approximation fails, the integrated likelihoods are all `NA`, which propagates
through the weight calculation. The final check:

```r
if (maci[2]/maci[1] > 20) { ... }
```

receives `NA / NA` = `NaN`, and `if (NaN)` throws the error.

## Affected location

| File | Approximate line | Context |
|------|-----------------|---------|
| `R/fun_Laplace.R` | MA credible-interval post-processing | `if (maci[2]/maci[1] > 20)` |

## Suggested fix

Guard against all-NA weights before the final post-processing:

```r
if (all(is.na(maci)) || is.na(maci[2]/maci[1])) {
  warning("All models failed the Laplace approximation; ",
          "no model-averaged BMD could be computed.")
  # Return NA BMD interval
  maci <- c(BMDL = NA_real_, BMD = NA_real_, BMDU = NA_real_)
} else if (maci[2]/maci[1] > 20) {
  # ... existing truncation logic ...
}
```

## Trigger conditions

The crash requires **all** models to fail, which happens when:

- Data is very sparse (few litters per dose group)
- The Beta-Binomial likelihood surface is flat or multimodal
- `rstan::optimizing()` cannot find a positive-definite Hessian for any model

With ≥ 8–10 litters per dose group the problem typically disappears.

## Workaround

Use larger datasets (≥ 10 litters per dose). In vignette examples we increased
from 3 to 10 litters per dose to avoid triggering this.

## Related

- `bug_report_gamlss_try_error.md` — unguarded `try()` on `gamlss` BB fit
- Both bugs affect the `cluster = TRUE` quantal pathway

## Environment

- BMABMDR 0.1.20
- R 4.5.1
- Windows 11
