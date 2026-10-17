## =============================================================================
## Helper functions for the BMABMDR simulation study
##
##   "BMD10 / BMD20 / BMD50: separate model-averaging fits (one per BMR)
##    versus a single fit of the best model used for all BMRs at once"
##
## See simulation/README.md for the design and simulation/PROGRESS.md for the
## work log. All functions are self-contained and only rely on exported
## BMABMDR functions (PREP_DATA_N, PREP_DATA_LN, full.laplace_MA, get_models,
## DRM.*).
## =============================================================================

## -----------------------------------------------------------------------------
## 1. True dose-response curves
## -----------------------------------------------------------------------------

#' Mean (Normal errors) / median (Lognormal errors) of the true curve.
#' All curves have background `a`; `c` is the maximal fold change mu(inf)/a.
true_curve <- function(x, truth) {
  p <- truth$par
  switch(truth$model,
         # exponential: a(1 + (c-1)(1 - exp(-b x^d)))
         E4  = p$a * (1 + (p$c - 1) * (1 - exp(-p$b * x^p$d))),
         # Hill: a(1 + (c-1) x^d / (b^d + x^d))
         H4  = p$a * (1 + (p$c - 1) * x^p$d / (p$b^p$d + x^p$d)),
         # log-normal (probit in log-dose)
         LN4 = p$a * (1 + (p$c - 1) * pnorm(p$d * (log(x) - log(p$b)))),
         # linear without plateau: misspecified for all 4-parameter models
         LIN = p$a * (1 + p$s * x),
         stop("unknown true model ", truth$model))
}

#' True BMD for a relative change `q` (BMABMDR definition: mu(BMD) = mu(0)(1 +/- q)).
#' Returns Inf when the true curve never reaches the BMR.
true_bmd <- function(truth, q, upper = 1e6) {
  mu0 <- true_curve(0, truth)
  inc <- true_curve(upper, truth) > mu0
  target <- if (inc) mu0 * (1 + q) else mu0 * (1 - q)
  f <- function(x) true_curve(x, truth) - target
  if (sign(f(upper)) == sign(f(0))) return(Inf)
  uniroot(f, c(0, upper), tol = 1e-12)$root
}

## -----------------------------------------------------------------------------
## 2. Data generation
## -----------------------------------------------------------------------------

#' Simulate individual responses and return dose-group summary statistics
#' (dose, arithmetic mean, sd, n) in the column order used by PREP_DATA_N/LN.
#'
#' error = "normal":    y ~ N(mu(x), (cv * mu(0))^2)  (constant variance)
#' error = "lognormal": log y ~ N(log mu(x), cv^2)     (mu = median)
simulate_summary_data <- function(truth, doses, n, cv, error = c("normal", "lognormal")) {
  error <- match.arg(error)
  mu <- true_curve(doses, truth)
  res <- lapply(seq_along(doses), function(i) {
    y <- if (error == "normal") {
      rnorm(n, mu[i], cv * mu[1])
    } else {
      exp(rnorm(n, log(mu[i]), cv))
    }
    c(x = doses[i], y = mean(y), s = sd(y), n = n)
  })
  as.data.frame(do.call(rbind, res))
}

## -----------------------------------------------------------------------------
## 3. BMD for an arbitrary BMR from posterior draws of a fitted BMABMDR model
## -----------------------------------------------------------------------------

#' Name of the DRM function for a model (e.g. "E4_N", increasing -> DRM.E4_NI)
drm_fun <- function(model, increasing) {
  get(paste0("DRM.", model, ifelse(increasing, "I", "D")),
      envir = asNamespace("BMABMDR"))
}

#' For each posterior draw (rows of `pars`, columns p1..p4 as returned in
#' `fit$parsN[[i]]` / `fit$parsLN[[i]]`) compute the BMD for BMR `q_new`, by
#' inverting the fitted curve that was parameterised with BMR `q_fit`.
#'
#' Doses are on the BMABMDR rescaled scale (dose / max dose).
#' Normal models: curve is the mean -> target mu(0)(1 +/- q_new).
#' Lognormal models: curve is log-median -> target log mu(0) + log(1 +/- q_new).
#' Draws for which the fitted asymptote does not reach the BMR get Inf;
#' draws for which the curve cannot be evaluated get NA.
bmd_from_draws <- function(pars, model, increasing, q_fit, q_new, shift = 0,
                           grid = c(0, exp(seq(log(1e-6), log(1e4), length.out = 400)))) {
  fun <- drm_fun(model, increasing)
  lognormal <- grepl("_LN$", model)
  P <- as.matrix(pars[, c("p1", "p2", "p3", "p4")])
  out <- rep(NA_real_, nrow(P))
  for (i in seq_len(nrow(P))) {
    par <- P[i, ]
    if (any(!is.finite(par))) next
    # asymptote check: p3 is the fold change mu(inf)/mu(0) for all models
    if (increasing && par[3] <= 1 + q_new) { out[i] <- Inf; next }
    if (!increasing && par[3] >= 1 - q_new) { out[i] <- Inf; next }
    fx <- suppressWarnings(fun(par, grid, q_fit, shift))
    if (!is.finite(fx[1])) next
    target <- if (lognormal) {
      fx[1] + log(if (increasing) 1 + q_new else 1 - q_new)
    } else {
      fx[1] * (if (increasing) 1 + q_new else 1 - q_new)
    }
    g <- if (increasing) fx - target else target - fx
    hit <- which(is.finite(g) & g >= 0)
    if (length(hit) == 0) { out[i] <- Inf; next }
    j <- hit[1]
    if (j == 1) { out[i] <- 0; next }
    h <- function(x) {
      v <- suppressWarnings(fun(par, x, q_fit, shift))
      if (increasing) v - target else target - v
    }
    r <- try(uniroot(h, c(grid[j - 1], grid[j]), tol = 1e-10)$root, silent = TRUE)
    if (!inherits(r, "try-error")) out[i] <- r
  }
  out
}

#' Quantiles that are robust to Inf values (non-estimable draws).
q_inf <- function(x, probs) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  unname(quantile(x, probs, type = 1))
}

## -----------------------------------------------------------------------------
## 4. Fitting wrappers
## -----------------------------------------------------------------------------

#' Prepare data and run full.laplace_MA for a single BMR, silencing the console
#' output of the package. Returns the fit (or a try-error) and run time.
fit_ma <- function(sumdata, q, prior.weights = rep(1, 16), ndraws = 30000,
                   seed = 123, pvec = c(0.05, 0.5, 0.95)) {
  t0 <- proc.time()[["elapsed"]]
  warns <- character()
  fit <- withCallingHandlers(
    try({
      utils::capture.output({
        dN  <- suppressMessages(PREP_DATA_N(sumdata, sumstats = TRUE, sd = TRUE, q = q))
        dLN <- suppressMessages(PREP_DATA_LN(sumdata, sumstats = TRUE, sd = TRUE, q = q))
        f <- suppressMessages(full.laplace_MA(dN, dLN, prior.weights = prior.weights,
                                              ndraws = ndraws, seed = seed, pvec = pvec))
      })
      f$is_increasing <- dN$data$is_increasing == 1
      f
    }, silent = TRUE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(fit = fit, time = proc.time()[["elapsed"]] - t0, warnings = unique(warns))
}

fit_ok <- function(fit) {
  !inherits(fit, "try-error") && !is.null(fit$weights) && all(is.finite(fit$MA))
}

#' Draws (data.frame with p1..p4, BMD) of model number `i` (get_models order).
model_draws <- function(fit, i) {
  if (i <= 8) fit$parsN[[i]] else fit$parsLN[[i - 8]]
}

## -----------------------------------------------------------------------------
## 5. Estimators
## -----------------------------------------------------------------------------
## All estimators return a data.frame with one row per BMR:
##   q, BMDL, BMD, BMDU (original dose scale), p_inf (proportion of
##   non-estimable posterior draws), best_model

#' (A) MA_sep: model averaging, one separate fit per BMR (package default use)
#' (B) Best_sep: best-weight model of each separate fit, its own BMD posterior
est_separate <- function(fits, qs, pvec) {
  ma <- best <- vector("list", length(qs))
  mods <- unname(BMABMDR::get_models("continuous"))
  for (k in seq_along(qs)) {
    f <- fits[[k]]$fit
    if (!fit_ok(f)) {
      ma[[k]] <- best[[k]] <- data.frame(q = qs[k], BMDL = NA, BMD = NA, BMDU = NA,
                                         p_inf = NA, best_model = NA)
      next
    }
    ma[[k]] <- data.frame(q = qs[k], BMDL = f$MA[1], BMD = f$MA[2], BMDU = f$MA[3],
                          p_inf = 0, best_model = NA)
    ib <- which.max(f$weights)
    dr <- model_draws(f, ib)
    bq <- q_inf(dr$BMD * f$max.dose, pvec)
    best[[k]] <- data.frame(q = qs[k], BMDL = bq[1], BMD = bq[2], BMDU = bq[3],
                            p_inf = 0, best_model = mods[ib])
  }
  list(MA_sep = do.call(rbind, ma), Best_sep = do.call(rbind, best))
}

#' (C) Best_joint: ONE fit at q_ref, best-weight model, all BMRs from its draws
#' (D) MA_joint: ONE fit at q_ref, all BMRs from the weight-mixture of all
#'     models' draws (model averaging of the derived BMDs)
est_joint <- function(fit_ref, q_ref, qs, pvec, n_draws = 2000, min_weight = 1e-4) {
  mods <- unname(BMABMDR::get_models("continuous"))
  empty <- data.frame(q = qs, BMDL = NA, BMD = NA, BMDU = NA, p_inf = NA, best_model = NA)
  if (!fit_ok(fit_ref)) return(list(Best_joint = empty, MA_joint = empty))
  inc <- fit_ref$is_increasing
  w <- fit_ref$weights
  w[!is.finite(w)] <- 0
  ib <- which.max(w)
  use <- which(w > min_weight)
  w_use <- w[use] / sum(w[use])
  counts <- as.vector(rmultinom(1, n_draws, w_use))
  # sub-sample draws per model once, so all BMRs use the same posterior draws
  sub <- lapply(seq_along(use), function(j) {
    dr <- model_draws(fit_ref, use[j])
    nn <- if (use[j] == ib) max(counts[j], n_draws) else counts[j]
    dr[sample.int(nrow(dr), nn, replace = TRUE), , drop = FALSE]
  })
  best_rows <- ma_rows <- vector("list", length(qs))
  for (k in seq_along(qs)) {
    bm <- lapply(seq_along(use), function(j) {
      if (nrow(sub[[j]]) == 0) return(numeric())
      bmd_from_draws(sub[[j]], mods[use[j]], inc, q_ref, qs[k], fit_ref$shift) *
        fit_ref$max.dose
    })
    jb <- which(use == ib)
    b_all <- bm[[jb]]
    # Best: all draws of the best model; MA: counts[j] draws of every model
    m_all <- unlist(lapply(seq_along(use), function(j) {
      if (use[j] == ib) head(bm[[j]], counts[j]) else bm[[j]]
    }))
    bq <- q_inf(b_all, pvec); mq <- q_inf(m_all, pvec)
    best_rows[[k]] <- data.frame(q = qs[k], BMDL = bq[1], BMD = bq[2], BMDU = bq[3],
                                 p_inf = mean(is.infinite(b_all)), best_model = mods[ib])
    ma_rows[[k]] <- data.frame(q = qs[k], BMDL = mq[1], BMD = mq[2], BMDU = mq[3],
                               p_inf = mean(is.infinite(m_all)), best_model = NA)
  }
  list(Best_joint = do.call(rbind, best_rows), MA_joint = do.call(rbind, ma_rows))
}

## -----------------------------------------------------------------------------
## 6. One replicate
## -----------------------------------------------------------------------------

#' Simulate one data set for scenario `sc` and apply all four estimators.
run_one_replicate <- function(sc, rep_id, qs, q_ref, pvec, ndraws, n_draws_joint,
                              prior.weights, seed) {
  set.seed(seed)
  dat <- simulate_summary_data(sc$truth, sc$doses, sc$n, sc$cv, sc$error)
  # separate fits (one per BMR); the fit at q_ref is re-used by the joint arms
  fits <- lapply(qs, function(q) fit_ma(dat, q, prior.weights, ndraws, seed, pvec))
  sep <- est_separate(fits, qs, pvec)
  t0 <- proc.time()[["elapsed"]]
  jnt <- est_joint(fits[[match(q_ref, qs)]]$fit, q_ref, qs, pvec, n_draws_joint)
  t_post <- proc.time()[["elapsed"]] - t0
  tb <- sapply(qs, function(q) true_bmd(sc$truth, q))
  times <- sapply(fits, `[[`, "time")
  t_ref <- times[match(q_ref, qs)]
  res <- do.call(rbind, lapply(names(c(sep, jnt)), function(m) {
    d <- c(sep, jnt)[[m]]
    d$method <- m
    d
  }))
  res$true_BMD <- tb[match(res$q, qs)]
  res$scenario <- sc$id
  res$rep <- rep_id
  res$fit_time <- ifelse(grepl("_sep$", res$method), sum(times), t_ref + t_post)
  attr(res, "data") <- dat
  attr(res, "warnings") <- lapply(fits, `[[`, "warnings")
  res
}

## -----------------------------------------------------------------------------
## 7. Performance measures
## -----------------------------------------------------------------------------

#' Per-replicate derived measures (added as columns).
add_measures <- function(res) {
  res$ok <- is.finite(res$BMD) & is.finite(res$BMDL)
  res$log_ratio <- log(res$BMD / res$true_BMD)            # log bias of median
  res$cover <- res$BMDL <= res$true_BMD & res$true_BMD <= res$BMDU
  res$bmdl_below <- res$BMDL <= res$true_BMD              # BMDL protective?
  res$log10_width <- log10(res$BMDU / res$BMDL)           # precision
  res
}

#' Within a data set and method: are BMD10 <= BMD20 <= BMD50 (and same for BMDL)?
ordering_check <- function(res) {
  sp <- split(res, list(res$scenario, res$rep, res$method), drop = TRUE)
  do.call(rbind, lapply(sp, function(d) {
    d <- d[order(d$q), ]
    data.frame(scenario = d$scenario[1], rep = d$rep[1], method = d$method[1],
               bmd_ordered = if (anyNA(d$BMD)) NA else !is.unsorted(d$BMD),
               bmdl_ordered = if (anyNA(d$BMDL)) NA else !is.unsorted(d$BMDL))
  }))
}

#' Aggregate over replicates, with Monte Carlo standard errors for coverage.
summarise_performance <- function(res) {
  res <- add_measures(res)
  grp <- split(res, list(res$scenario, res$method, res$q), drop = TRUE)
  out <- do.call(rbind, lapply(grp, function(d) {
    ok <- d$ok
    nr <- nrow(d)
    cov <- mean(d$cover[ok])
    data.frame(
      scenario = d$scenario[1], method = d$method[1], q = d$q[1],
      true_BMD = d$true_BMD[1], n_rep = nr,
      fail_rate = mean(!ok),
      median_log_ratio = median(d$log_ratio[ok]),
      mean_log_ratio = mean(d$log_ratio[ok]),
      rmse_log = sqrt(mean(d$log_ratio[ok]^2)),
      coverage = cov, coverage_mcse = sqrt(cov * (1 - cov) / sum(ok)),
      bmdl_below_true = mean(d$bmdl_below[ok]),
      median_log10_width = median(d$log10_width[ok]),
      mean_p_inf = mean(d$p_inf, na.rm = TRUE),
      mean_fit_time = mean(d$fit_time, na.rm = TRUE)
    )
  }))
  rownames(out) <- NULL
  out[order(out$scenario, out$q, out$method), ]
}
