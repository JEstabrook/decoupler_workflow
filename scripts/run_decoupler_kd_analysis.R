library(decoupleR)
library(decoupleRBench)
library(dplyr)
library(tibble)

args = commandArgs(trailingOnly=TRUE)

expr_fname = args[1]
meta_fname = args[2]
netw_fname = args[3]
knock_down = args[4]
celltype = args[5]
results_out = args[6]
temp_expr = args[7]
temp_meta = args[8]
temp_netw = args[9]
pert_time = args[10]

expr <- readRDS(expr_fname)
meta <- readRDS(meta_fname)
network <- readRDS(netw_fname)

regulators <- c("ERBB2","MAP2K1","PIK3CG","EGFR","BCR","ERBB4","CDK1","PIK3CA","MET","CDK9","BRAF","SRC","PLK1","CDK6","AKT1","PTK2","PDGFRB","YES1","MTOR","AURKB","PRKCB","PDGFRA","PIM1","TEK","CAMK2G","AURKA","MAPK12","TYRO3","RET","PARP1","BCL2L2","HDAC6","HDAC1","BRD3","HDAC4","PTGS2","DNMT1","MAPK14","TGFB1","ALK","IKBKB","TOP2A","MAPK7","GSK3B","MAP2K5","AKT2","KIT","CSF1R","MAPK3","CDK7","FGFR4","ATM","SIRT1","HDAC3","EHMT2","HDAC2","MAOB","ATR","PRKDC","CHEK2","RPS6KA3","RAF1","TTK","MAPK11","MAPK8","JAK2","AXL","MDM2","JAK3")

sub_meta <- meta[(meta$pert_id == eval(knock_down)) & (meta$cell_id == eval(celltype)) & (meta$target %in% regulators) & (meta$pert_time == eval(pert_time)),]
control_meta <- meta[(meta$pert_id == 'DMSO') & (meta$cell_id == eval(celltype)) & (meta$pert_time == eval(pert_time)),]
sub_control <- control_meta[control_meta$det_plate %in% sub_meta$det_plate,]
sub_control$target <- unique(sub_meta$target)
joined_meta <- bind_rows(sub_meta,sub_control)

sub_expr <- expr[,joined_meta$id]
sub_netw <- network[network$tf %in% regulators,]

saveRDS(joined_meta,file=temp_meta)
saveRDS(sub_expr,file=temp_expr)
saveRDS(sub_netw,file=temp_netw)

seed <- 1
nproc = 4

stats_list = list(c('aucell','wmean','wsum','ulm','viper','gsva','ora','fgsea','udt','mdt','enricher'))
opts_list <- list(list(
  udt = list(sparse=FALSE, center=FALSE, na.rm=FALSE, min_n = 20, seed=seed),
  mdt = list(sparse=FALSE, center=FALSE, na.rm=FALSE, trees = 10, min_n = 20,
             nproc = nproc, seed=seed),
  aucell = list(nproc=nproc, seed=seed),
  wmean = list(times=100, sparse=TRUE, randomize_type = "rows", seed=seed),
  wsum = list(times=100, sparse=TRUE, randomize_type = "rows", seed=seed),
  ulm = list(sparse=FALSE, center=FALSE, na.rm=FALSE),
  viper = list(verbose=FALSE, pleiotropy=T, eset.filter=F),
  gsva = list(verbose = FALSE, method = "gsva"),
  ora = list(n_up=300, n_bottom=300, n_background=20000, with_ties = TRUE),
  fgsea = list(times=100, nproc=nproc, seed=seed),
  enricher = list(scaler_type="robust", minsize=5)
))

# Design
design <- tibble(
  set_name = 'pathway_commons', # name of the set resource
  bench_name = "knocktf", # name of the benchmark data
  stats_list = stats_list,
  opts_list = opts_list,
  bexpr_loc = temp_expr, # benchmark data location
  bmeta_loc = temp_meta, # metadata location
  source_loc = temp_netw, # set source location
  source_col = "tf", # source name of the gene set source
  target_col = "target", # target name of the set source
  filter_col = "confidence", # column by which we wish to filter
  filter_crit = list(c('A','B','C')) # criteria by which we wish to filter
)

# Run benchmark
result <- run_benchmark(
  .design = design, # provide input tibble
  .minsize = 5, # filter gene sets with size < 10
  .form = TRUE, # format the benchmark results
  .perform = TRUE, # evaluate benchmarking performance
  .silent = FALSE, # silently run the pipeline
  .downsample_pr = TRUE, # downsample TNs for precision-recall curve
  .downsample_roc = TRUE, # downsample TNs for ROC
  .downsample_times = 100, # downsampling iterations
  .url_bool = FALSE # whether to load from url
)

# Save result
saveRDS(result,file=results_out) 
