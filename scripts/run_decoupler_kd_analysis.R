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
enrich_slot = args[10]

expr <- readRDS(expr_fname)
meta <- readRDS(meta_fname)
network <- readRDS(netw_fname)

regulators <- c("STAT5A","BACH1", "PCBP1", "NR4A1", "FOXM1", "RCOR1", "GATA1","NFE2L1","MAZ", "HMGA1", "HDAC8", "STAT2", "E2F4","NFE2L2","POLR2G","SMARCA4","AGO1","KAT2B", "SMAD5", "NR2F2", "MITF","LMNA","TRIM28","E2F6","ATF3","MAFG","USF2","CEBPZ","CITED2","CTBP1", "NRF1","TBL1XR1","AGO2","GTF2F1","RELA","MAX", "NFYB","TFDP1", "JUND","STAT6", "STAT1", "SP1","NFATC1","SRF", "SNW1","ERF", "GATA2", "NFYA","BRCA1","HSF1","SRSF3")

sub_meta <- meta[(meta$knockdown_method == eval(knock_down)) & (meta$cell == eval(celltype)) & (meta$target %in% regulators),]
sub_expr <- expr[,sub_meta$id]
sub_netw <- network[network$tf %in% regulators,]

saveRDS(sub_meta,file=temp_meta)
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
  enricher = list(scaler_type="robust", minsize=5,enr_type=enrich_slot)
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
