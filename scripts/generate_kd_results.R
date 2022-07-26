library(dplyr)

args <- commandArgs(trailingOnly=TRUE)

res_handle <- args[1]
out_file <- args[2]
summary_file <- args[3]

results <- readRDS(res_handle)
joined_results <- bind_rows(results@bench_res$activity)
auc_sum <- bind_rows(results@summary$summary_table)

write.table(joined_results,out_file,sep='\t',quote=F,row.names=F)
write.table(auc_sum, summary_file,sep='\t',quote=F,row.names=F)
