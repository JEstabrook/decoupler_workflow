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

regulators<-c("ACTL6A","AHR","AR","ARNT","ATF2","ATF3","ATG5","ATM","ATR","BACH1","BRCA1","BRD4","CDX2","CHEK1","CREB1","EGR3","ELK1","EP300","EPAS1","EREG","ESR1","ESR2","ESRRA","ETS1","ETS2","ETV4","FOXA1","FOXM1","FOXO1","FOXO3","FOXP1","GATA2","GATA3","GATA6","GLI2","GTF3A","HES1","HIF1A","HMGA1","HNF1B","HOXA5","HSF1","HSF2","IRF2","IRF7","JUN","KLF4","LEF1","MAF","MITF","MSX1","MTF1","MYB","MYBL2","MYC","MYCN","NANOG","NFATC3","NFE2L2","NFKB1","NFKB2","NFYA","NFYC","NR2F2","NR3C1","PAX2","PAX3","PGR","PITX2","POSTN","POU5F1","PROX1","PTEN","RELA","RELB","REST","RUNX1","RUNX2","SF1","SMAD2","SMAD3","SMAD4","SNAI1","SOX2","SOX9","SP3","STAT1","STAT3","STAT5A","STAT6","TCF4","TCF7L1","TFAP2C","TFAP4","THRA","THRB","TP53","TP63","TWIST1","TYK2","XBP1","YY1","ZEB1","ZIC2")

sub_meta_ <- meta[(meta$knockdown_method == eval(knock_down)) & (meta$cell == eval(celltype)) & (meta$target %in% regulators),]

n_ <- dim(sub_meta_)[[1]]
control_sub_meta <- meta[(meta$knockdown_method == "Control")& (meta$cell == eval(celltype)),] 
cn_ <- dim(control_sub_meta)[[1]]
if(cn_ == 0){
    control_sample_meta <- meta[(meta$knockdown_method == "Control"),]
    sample_control <- control_sample_meta[sample(nrow(control_sample_meta),n_),]
} else {
    sample_control <- control_sub_meta}

#} else {
#    if(cn_ < n_){
#        dif_n <- n_ - cn_
#        control_sample_meta <- meta[(meta$knockdown_method == "Control"),]
#        sub_sample_control <- control_sample_meta[sample(nrow(control_sample_meta),dif_n),]
#        sample_control <- rbind(control_sub_meta,sub_sample_control)
#    } else {
#        sample_control <- control_sub_meta}
#}

n_targets <- length(unique(sub_meta_$target))
if(n_targets == 1){
    target <- unique(sub_meta_$target)
    length.n <- dim(sample_control)[1]
    sub_targs <- regulators[regulators != target]
    control_targets <- rep(sub_targs,length.out=length.n)
    sample_control$target <- control_targets
} else {
    length.n <- dim(sample_control)[1]
    targets <- unique(sub_meta_$target)
    sub_targs <- regulators[!(regulators %in% targets)]
    control_targets <- rep(sub_targs,length.out=length.n)
    sample_control$target <- control_targets
}
sub_meta <- rbind(sub_meta_,sample_control)

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
