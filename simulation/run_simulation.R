## =============================================================================
## Driver: BMABMDR simulation study
##   BMD10/20/50 by separate model-averaging fits vs. one best-model fit
##
## Usage (from the repository root):
##   Rscript simulation/run_simulation.R                 # pilot (default)
##   SIM_MODE=full SIM_CORES=8 Rscript simulation/run_simulation.R
##   SIM_MODE=full SIM_NREP=100 Rscript simulation/run_simulation.R
##
## Environment variables (all optional):
##   SIM_MODE     "pilot" (2 scenarios x few reps) or "full" (all scenarios)
##   SIM_NREP     replicates per scenario (pilot: 3, full: 200)
##   SIM_CORES    parallel workers (forked; default: detectCores() - 1)
##   SIM_NDRAWS   Laplace draws per model in full.laplace_MA (default 10000)
##   SIM_QREF     reference BMR for the single-fit arms (default 0.1)
##   SIM_OUT      output directory (default simulation/results/<mode>)
##
## Each replicate is saved as its own .rds file in <out>/raw/; replicates that
## already exist are skipped, so an interrupted run is resumed simply by
## re-running the same command (see simulation/PROGRESS.md).
## =============================================================================

suppressPackageStartupMessages({
  library(BMABMDR)
  library(parallel)
})

args_dir <- if (file.exists("simulation/R/sim_functions.R")) "simulation" else "."
source(file.path(args_dir, "R", "sim_functions.R"))

env <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (v == "") default else v
}
mode    <- env("SIM_MODE", "pilot")
nrep    <- as.integer(env("SIM_NREP", if (mode == "full") 200 else 3))
cores   <- as.integer(env("SIM_CORES", max(1, parallel::detectCores() - 1)))
ndraws  <- as.integer(env("SIM_NDRAWS", 10000))
q_ref   <- as.numeric(env("SIM_QREF", 0.1))
out_dir <- env("SIM_OUT", file.path(args_dir, "results", mode))

qs   <- c(0.1, 0.2, 0.5)
pvec <- c(0.05, 0.5, 0.95)            # 90% credible interval, as in BMABMDR
n_draws_joint <- 2000                 # posterior draws used for the inversion
prior.weights <- rep(1, 16)           # all 16 continuous models
stopifnot(q_ref %in% qs)

## -----------------------------------------------------------------------------
## Scenario grid
## -----------------------------------------------------------------------------
doses <- c(0, 6.25, 12.5, 25, 50, 100)

truths <- list(
  E4_inc  = list(model = "E4",  par = list(a = 10, b = 0.0144, c = 3,   d = 1)),
  H4_inc  = list(model = "H4",  par = list(a = 10, b = 30,     c = 3,   d = 2)),
  LN4_inc = list(model = "LN4", par = list(a = 10, b = 40,     c = 3,   d = 1.5)),
  LIN_inc = list(model = "LIN", par = list(a = 10, s = 0.02)),
  H4_dec  = list(model = "H4",  par = list(a = 10, b = 30,     c = 0.3, d = 2))
)

grid <- expand.grid(truth = names(truths), error = c("normal", "lognormal"),
                    cv = c(0.10, 0.20), n = 10, stringsAsFactors = FALSE)
if (mode == "pilot") {
  grid <- grid[grid$truth %in% c("H4_inc", "H4_dec") & grid$error == "normal" &
                 grid$cv == 0.10, ]
}
grid$id <- sprintf("%s_%s_cv%02d_n%d", grid$truth, grid$error, round(100 * grid$cv), grid$n)
grid$scen_no <- seq_len(nrow(grid))

scenarios <- lapply(seq_len(nrow(grid)), function(i) {
  list(id = grid$id[i], truth = truths[[grid$truth[i]]], doses = doses,
       n = grid$n[i], cv = grid$cv[i], error = grid$error[i], scen_no = grid$scen_no[i])
})

dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)

## true BMDs for the record
truth_tab <- do.call(rbind, lapply(scenarios, function(sc) {
  data.frame(scenario = sc$id, q = qs, true_BMD = sapply(qs, function(q) true_bmd(sc$truth, q)))
}))
write.csv(truth_tab, file.path(out_dir, "true_bmds.csv"), row.names = FALSE)
saveRDS(list(grid = grid, truths = truths, doses = doses, qs = qs, q_ref = q_ref,
             pvec = pvec, ndraws = ndraws, n_draws_joint = n_draws_joint,
             prior.weights = prior.weights, nrep = nrep),
        file.path(out_dir, "settings.rds"))

## -----------------------------------------------------------------------------
## Run (one task = one scenario x replicate; checkpointed per task)
## -----------------------------------------------------------------------------
tasks <- expand.grid(s = seq_along(scenarios), r = seq_len(nrep))
tasks$file <- file.path(out_dir, "raw",
                        sprintf("%s_rep%03d.rds", sapply(tasks$s, function(i) scenarios[[i]]$id), tasks$r))
todo <- which(!file.exists(tasks$file))
message(sprintf("[%s] mode=%s: %d scenarios x %d reps = %d tasks, %d to do, %d cores",
                format(Sys.time()), mode, length(scenarios), nrep, nrow(tasks), length(todo), cores))

run_task <- function(i) {
  sc <- scenarios[[tasks$s[i]]]
  r  <- tasks$r[i]
  seed <- 20261017L + 10000L * sc$scen_no + r
  res <- try(run_one_replicate(sc, r, qs, q_ref, pvec, ndraws, n_draws_joint,
                               prior.weights, seed), silent = TRUE)
  if (inherits(res, "try-error")) {
    cat(sprintf("task %d (%s rep %d) failed: %s\n", i, sc$id, r, res), file = stderr())
    return(FALSE)
  }
  saveRDS(res, tasks$file[i])
  TRUE
}

t0 <- Sys.time()
ok <- if (cores > 1) {
  unlist(mclapply(todo, run_task, mc.cores = cores, mc.preschedule = FALSE))
} else {
  vapply(todo, run_task, logical(1))
}
message(sprintf("[%s] finished %d/%d tasks in %.1f min", format(Sys.time()), sum(ok),
                length(todo), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## -----------------------------------------------------------------------------
## Collect
## -----------------------------------------------------------------------------
files <- list.files(file.path(out_dir, "raw"), pattern = "\\.rds$", full.names = TRUE)
all_res <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f)
  attr(x, "data") <- NULL; attr(x, "warnings") <- NULL
  x
}))
saveRDS(all_res, file.path(out_dir, "all_results.rds"))
message("collected ", length(files), " replicates -> ", file.path(out_dir, "all_results.rds"))
message("next: Rscript simulation/analyze_results.R ", out_dir)
