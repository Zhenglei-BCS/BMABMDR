## =============================================================================
## Analysis of the BMABMDR simulation study results
##
## Usage (from the repository root):
##   Rscript simulation/analyze_results.R [results_dir]
## default results_dir: simulation/results/pilot
##
## Writes to <results_dir>:
##   performance_summary.csv   per scenario x BMR x method
##   ordering_summary.csv      % of data sets with BMD10<=BMD20<=BMD50 per method
##   overall_summary.csv       averaged over scenarios, per BMR x method
##   best_model_agreement.csv  how often the best model differs between the separate fits
##   fig_*.png                 figures
## =============================================================================

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
base <- if (file.exists("simulation/R/sim_functions.R")) "simulation" else "."
out_dir <- if (length(args) > 0) args[1] else file.path(base, "results", "pilot")
source(file.path(base, "R", "sim_functions.R"))

res <- readRDS(file.path(out_dir, "all_results.rds"))
res$method <- factor(res$method, levels = c("MA_sep", "MA_joint", "Best_sep", "Best_joint"))
res$BMR <- factor(paste0("BMD", round(100 * res$q)), levels = c("BMD10", "BMD20", "BMD50"))

## performance per scenario
perf <- summarise_performance(res)
write.csv(perf, file.path(out_dir, "performance_summary.csv"), row.names = FALSE)

## overall (scenario-averaged)
overall <- aggregate(cbind(fail_rate, median_log_ratio, rmse_log, coverage,
                           bmdl_below_true, median_log10_width, mean_p_inf,
                           mean_fit_time) ~ method + q,
                     data = perf, FUN = mean)
overall <- overall[order(overall$q, overall$method), ]
write.csv(overall, file.path(out_dir, "overall_summary.csv"), row.names = FALSE)

## ordering consistency
ord <- ordering_check(res)
ord_sum <- aggregate(cbind(bmd_ordered, bmdl_ordered) ~ method + scenario, data = ord,
                     FUN = function(x) mean(x), na.action = na.omit)
write.csv(ord_sum, file.path(out_dir, "ordering_summary.csv"), row.names = FALSE)

## does the best model change with the BMR in the separate fits?
bs <- res[res$method == "Best_sep", ]
agree <- do.call(rbind, lapply(split(bs, list(bs$scenario, bs$rep), drop = TRUE), function(d) {
  data.frame(scenario = d$scenario[1], rep = d$rep[1],
             same_best_all_q = length(unique(na.omit(d$best_model))) == 1)
}))
agree_sum <- aggregate(same_best_all_q ~ scenario, data = agree, FUN = mean)
write.csv(agree_sum, file.path(out_dir, "best_model_agreement.csv"), row.names = FALSE)

print(overall, digits = 3, row.names = FALSE)
cat("\nOrdering consistency (share of data sets with BMD10 <= BMD20 <= BMD50):\n")
print(aggregate(cbind(bmd_ordered, bmdl_ordered) ~ method, data = ord, FUN = mean), digits = 3)
cat("\nShare of data sets where the best model is the same at all BMRs (separate fits):\n")
print(agree_sum, digits = 3)

## -----------------------------------------------------------------------------
## Figures
## -----------------------------------------------------------------------------
m <- add_measures(res)
m <- m[m$ok, ]

p1 <- ggplot(m, aes(x = method, y = log_ratio, fill = method)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.size = 0.5) +
  facet_grid(BMR ~ scenario, scales = "free_y") +
  labs(y = "log(BMD estimate / true BMD)", x = NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                     legend.position = "none", strip.text.x = element_text(size = 6))
ggsave(file.path(out_dir, "fig_log_ratio.png"), p1, width = 14, height = 7, dpi = 120)

perf$BMR <- factor(paste0("BMD", round(100 * perf$q)), levels = c("BMD10", "BMD20", "BMD50"))
perf$method <- factor(perf$method, levels = levels(res$method))
p2 <- ggplot(perf, aes(x = BMR, y = coverage, colour = method, group = method)) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  geom_point(position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = coverage - 2 * coverage_mcse, ymax = coverage + 2 * coverage_mcse),
                width = 0, position = position_dodge(0.4)) +
  facet_wrap(~ scenario) + ylim(0, 1) +
  labs(y = "Coverage of 90% credible interval (BMDL, BMDU)", x = NULL) + theme_bw()
ggsave(file.path(out_dir, "fig_coverage.png"), p2, width = 12, height = 8, dpi = 120)

p3 <- ggplot(perf, aes(x = BMR, y = median_log10_width, colour = method, group = method)) +
  geom_point(position = position_dodge(0.4)) + geom_line(position = position_dodge(0.4)) +
  facet_wrap(~ scenario) +
  labs(y = "median log10(BMDU / BMDL)", x = NULL) + theme_bw()
ggsave(file.path(out_dir, "fig_width.png"), p3, width = 12, height = 8, dpi = 120)

cat("\nOutputs written to", out_dir, "\n")
