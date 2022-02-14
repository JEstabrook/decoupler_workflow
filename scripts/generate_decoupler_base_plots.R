library(decoupleR)
library(decoupleRBench)

args = commandArgs(trailingOnly=TRUE)

results_f = args[1]
roc_plot = args[2]
pr_plot = args[3]
auroc_heat  = args[4]
pr_heat = args[5]

results <- readRDS(results_f)

pdf(roc_plot)
results@summary$roc_plot
dev.off()

pdf(pr_plot)
results@summary$pr_plot
dev.off()

pdf(auroc_heat)
results@summary$auroc_heat
dev.off()

pdf(pr_heat)
results@summary$pr_heat
dev.off()

