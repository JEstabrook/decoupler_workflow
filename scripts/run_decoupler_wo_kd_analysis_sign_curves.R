library(decoupleR)
library(decoupleRBench)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyselect)
library(purrr)
library(yardstick)
library(stringr)
library(ggplot2)
library(pheatmap)

args = commandArgs(trailingOnly=TRUE)

expr_fname = args[1]
meta_fname = args[2]
netw_fname = args[3]
celltype = args[4]
results_out = args[5]
temp_expr = args[6]
temp_meta = args[7]
temp_netw = args[8]
enrich_slot = args[9]

expr <- readRDS(expr_fname)
meta <- readRDS(meta_fname)
network <- readRDS(netw_fname)

regulators<-c("ACTL6A","AHR","AR","ARNT","ATF2","ATF3","ATG5","ATM","ATR","BACH1","BRCA1","BRD4","CDX2","CHEK1","CREB1","EGR3","ELK1","EP300","EPAS1","EREG","ESR1","ESR2","ESRRA","ETS1","ETS2","ETV4","FOXA1","FOXM1","FOXO1","FOXO3","FOXP1","GATA2","GATA3","GATA6","GLI2","GTF3A","HES1","HIF1A","HMGA1","HNF1B","HOXA5","HSF1","HSF2","IRF2","IRF7","JUN","KLF4","LEF1","MAF","MITF","MSX1","MTF1","MYB","MYBL2","MYC","MYCN","NANOG","NFATC3","NFE2L2","NFKB1","NFKB2","NFYA","NFYC","NR2F2","NR3C1","PAX2","PAX3","PGR","PITX2","POSTN","POU5F1","PROX1","PTEN","RELA","RELB","REST","RUNX1","RUNX2","SF1","SMAD2","SMAD3","SMAD4","SNAI1","SOX2","SOX9","SP3","STAT1","STAT3","STAT5A","STAT6","TCF4","TCF7L1","TFAP2C","TFAP4","THRA","THRB","TP53","TP63","TWIST1","TYK2","XBP1","YY1","ZEB1","ZIC2")

sub_meta_ <- meta[(meta$cell == eval(celltype)) & (meta$target %in% regulators),]
n_ <- dim(sub_meta_)[[1]]
control_sub_meta <- meta[(meta$knockdown_method == "Control")& (meta$cell == eval(celltype)),]
cn_ <- dim(control_sub_meta)[[1]]
if(cn_ == 0){
    control_sample_meta <- meta[(meta$knockdown_method == "Control"),]
    sample_control <- control_sample_meta[sample(nrow(control_sample_meta),n_),]
} else {
    if(cn_ < n_){
        dif_n <- n_ - cn_
        control_sample_meta <- meta[(meta$knockdown_method == "Control"),]
        sub_sample_control <- control_sample_meta[sample(nrow(control_sample_meta),dif_n),]
        sample_control <- rbind(control_sub_meta,sub_sample_control)
    } else {
        sample_control <- control_sub_meta}
}

n_targets <- length(unique(sub_meta_$target))
if(n_targets == 1){
    sample_control$target <- unique(sub_meta_$target)
} else {
    length.n <- dim(sample_control)[1]
    targets <- unique(sub_meta_$target)
    control_targets <- rep(targets,length.out=length.n)
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


run_benchmark_mod <- function(.design,
                          .form = TRUE,
                          .perform = TRUE,
                          .minsize = 10,
                          .silent = TRUE,
                          .downsample_pr = FALSE,
                          .downsample_roc = FALSE,
                          .downsample_times = 100,
                          .url_bool = FALSE
){

  bench_env <- new.env()
  `%!in%` <- Negate(`%in%`)
  .design_formatted <- .design %>% decoupleRBench::format_design()
  .design_final <- .design_formatted[,names(.design_formatted)  %!in% c('extra_name','consensus_crit')]
  print(.design_final)
  res <- .design_final %>%
    dplyr::mutate(activity = pmap(.l=.,
                           .f=function(set_name, bench_name,
                                       stats_list, opts_list,
                                       bexpr_loc, bmeta_loc, source_loc,
                                       source_col, target_col, filter_col,
                                       filter_crit, noise_crit, weight_crit,
                                       .source_bln, .expr_bln, .meta_bln){

                             # Check prerequisites
                             if(!.expr_bln){
                               bench_env$gene_expression <- decoupleRBench::readRDS_helper(bexpr_loc, .url_bool) %>%
                                 as.matrix()
                               message(stringr::str_glue("Expression loaded"))
                             }
                             if(!.meta_bln){
                               bench_env$meta_data <- decoupleRBench::readRDS_helper(bmeta_loc, .url_bool)
                               message(stringr::str_glue("Meta loaded"))
                             }

                             # Filter set_source/network
                             bench_env$set_source <- decoupleRBench::readRDS_helper(source_loc, .url_bool)
                             ss_filtered <- filter_sets(bench_env$set_source, source_col,
                                                        filter_col, filter_crit,
                                                        .minsize, .silent)
                             message(stringr::str_glue("Network loaded"))

                             # Add noise
                             if (is.list(noise_crit)){
                               message(
                                 stringr::str_glue("Modify network"))
                               ss_filtered <- decoupleRBench::net_noise(
                                 network = ss_filtered,
                                 mode = noise_crit$mode,
                                 perc = noise_crit$perc,
                                 seed = noise_crit$seed,
                                 source = source_col,
                                 target = target_col
                               )
                               message(
                                 stringr::str_glue("{noise_crit$mode} {noise_crit$perc} noise"))
                             }

                             # Remove weight
                             if (is.list(weight_crit)){
                               message(
                                 stringr::str_glue("Unweight network"))
                               ss_filtered <- decoupleRBench::net_weight(
                                 network = ss_filtered,
                                 weight_rm = c(weight_crit$.mor, weight_crit$.likelihood)
                               )
                               message(
                                 stringr::str_glue("{c(weight_crit$.mor, weight_crit$.likelihood)} set to 1 "))
                             }

                             # Show Current Row/Run
                             if(!.silent){
                               .curr_row <- paste(set_name, bench_name,
                                                  paste0(unlist(filter_crit), collapse=""),
                                                  sep="_")
                               message(str_glue("Currently Running: {.curr_row}"))
                             }

                             # Match target genes between network and mat
                             targets <- rownames(bench_env$gene_expression)
                             msk <- ss_filtered[[target_col]] %in% targets
                             ss_filtered <- ss_filtered[msk,]
                             ss_filtered <- ss_filtered %>%
                               group_by_at(source_col) %>%
                               filter(n() >= .minsize)
                             # Obtain Activity with decouple and format
                             decoupleR::decouple(mat = bench_env$gene_expression, network = ss_filtered,
                                      .source = source_col, .target = tidyselect::all_of(target_col),
                                      statistics = stats_list, args = opts_list,
                                      include_time = TRUE)  %>%
                               dplyr::rename(id=.data$condition) %>%
                               ungroup() %>%
                               dplyr::left_join(bench_env$meta_data, by="id")  %>%
                               dplyr::group_by(.data$statistic) %>%
                               dplyr::group_split(.keep=TRUE) %>%
                               as.list()
                           })) %>% {
                             if(.form & !.perform) bench_format(., .silent)
                             else if(.form & .perform) bench_format(., .silent) %>%
                               mutate(roc = .data$activity %>%
                                        map(~calc_curve_mod(df=.x,
                                                        downsampling=.downsample_roc,
                                                        times=.downsample_times,
                                                        curve="ROC")),
                                      prc = .data$activity %>%
                                        map(~calc_curve_mod(df=.x,
                                                        downsampling=.downsample_pr,
                                                        times=.downsample_times,
                                                        curve="PR")))
                             else .
                           }

  if(.form & .perform){
    bench_result <- methods::new("BenchResult",
                                 bench_res=res,
                                 summary=res %>% get_bench_summary(),
                                 design=.design)
  }
  else{
    bench_result <- methods::new("BenchResult",
                                 bench_res=res,
                                 summary=list(NULL),
                                 design=.design)
  }
  return(bench_result)
}

filter_sets <- function(set_source,
                        source_col,
                        filter_col,
                        filter_crit,
                        .minsize,
                        .silent){
  n_duprows <- sum(duplicated(set_source))
  na_bool <- is.na(filter_col)

  gs_filtered <- set_source %>%
    {
      if(na_bool){distinct(.)}
      else if(!na_bool){
        filter(., .data[[filter_col]] %in% filter_crit) %>%
          distinct_at(vars(-.data[[filter_col]]), .keep_all = FALSE)
      }
    } %>%
    group_by(.data[[source_col]]) %>%
    add_count() %>%
    filter(n >= .minsize) %>%
    ungroup()

  if (n_duprows & !.silent){
    warning(str_glue("{n_duprows} rows were duplicated in the set resource! ",
                     "{sum(duplicated(gs_filtered))} duplicated rows ",
                     "remain after filtering."))
  }
  return(gs_filtered)
}

calc_curve_mod = function(df,
                      downsampling = FALSE,
                      times = 1000,
                      curve = "ROC",
                      seed = 420){
  set.seed(seed)

  if(curve=="PR"){
    res_col_1 <- "precision"
    res_col_2 <- "recall"
    curve_fun = yardstick::pr_curve
    auc_fun = yardstick::pr_auc
  }
  else if(curve=="ROC"){
    res_col_1 <- "sensitivity"
    res_col_2 <- "specificity"
    curve_fun = yardstick::roc_curve
    auc_fun = yardstick::roc_auc
  }

  df = df %>%
    prepare_for_roc_curve(., filter_tn = TRUE)


  if (sum(which(df$response == 0)) == nrow(df)){
    return(as_tibble(NULL))
  }

  cn = df %>% filter(.data$response == 0)
  cp = df %>% filter(.data$response == 1)

  feature_coverage = length(unique(df$source))

  if (downsampling == TRUE) {
    num_tp = nrow(cp)

    res = map_df(seq(from=1, to=times, by=1), function(i) {
      df_sub = sample_n(cn, num_tp, replace=TRUE) %>%
        bind_rows(cp)

      r_sub = df_sub %>%
        curve_fun(.data$response, .data$predictor)

      auc = df_sub %>%
        auc_fun(.data$response, .data$predictor) %>%
        pull(.data$.estimate)

      res_sub = tibble({{ res_col_1 }} := r_sub %>% pull(res_col_1),
                       {{ res_col_2 }} := r_sub %>% pull(res_col_2),
                       th = r_sub$.threshold,
                       auc = auc,
                       n = length(which(df$response == 1)),
                       cp = nrow(cp),
                       cn = nrow(cn),
                       coverage = feature_coverage) %>%
        mutate("run" = i)

    })
    # Get Average AUC
    res <- res %>% dplyr::rename("raw_auc" = auc)
    # auc is the mean of all iterations, raw_auc is the value per iteration
    res$auc <- sum(res$raw_auc)/length(res$raw_auc)
    res$cn <- nrow(cp)

  } else {
    r = df %>%
      curve_fun(.data$response, .data$predictor)
    auc = df %>%
      auc_fun(.data$response, .data$predictor)

    res = tibble({{ res_col_1 }} := r %>% pull(res_col_1),
                 {{ res_col_2 }} := r %>% pull(res_col_2),
                 th = r$.threshold,
                 auc = auc$.estimate,
                 n = length(which(df$response == 1)),
                 cp = nrow(cp),
                 cn = nrow(cn),
                 coverage = feature_coverage) %>%
      arrange(!!res_col_1, !!res_col_2)
  }

  return(res)
}


prepare_for_roc_curve = function(df, filter_tn = FALSE) {
  res = df %>%
    dplyr::mutate(response = case_when(.data$source == .data$target & .data$sign == -1 ~ 1,
                                       .data$source == .data$target & .data$sign == 1 ~ 0),
                  predictor = .data$score)#*sign)
  res$response = factor(res$response, levels = c(1, 0))

  if (filter_tn == TRUE) {
    z = intersect(res$source, res$target)
    res = res %>%
      filter(.data$source %in% z, .data$target %in% z)
  }
  res %>%
    select(.data$source, .data$id, .data$response, .data$predictor)
}


get_bench_summary <- function(.res_tibble) {
  # get roc results
  roc <- format_roc(.res_tibble, "roc")

  # get PR roc results
  pr <- format_roc(.res_tibble, "prc")

  # Plot ROC
  roc_plot <- ggplot(roc, aes(x = 1-specificity,
                              y = sensitivity,
                              colour = .data$run_key)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")

  # Plot PR ROC
  # Assign Expected Random PR baseline (C.Positives/P+N)
  pr_plot <- ggplot(pr, aes(x = recall, y = precision, colour = .data$run_key)) +
    geom_line() +
    geom_abline(pr %>% rename("random_baseline"=.data$name_lvl),
                mapping=aes(intercept = .data$cp/(.data$cn+.data$cp), slope = 0,
                            linetype = .data$random_baseline)) +
    xlab("Recall/Sensitivity") +
    ylab("Precision")

  # Extract AUROC
  auroc_tibble <- .res_tibble %>%
    unnest(roc) %>%
    select(.data$set_name, .data$bench_name, .data$filter_crit,
           .data$statistic, .data$auc) %>%
    distinct()

  # Plot AUROC
  auroc_plot <- auroc_tibble %>%
    unite("run_key", .data$set_name, .data$bench_name,
          .data$statistic, .data$filter_crit, remove = FALSE) %>%
    ggplot(., aes(x = reorder(.data$run_key, .data$auc),
                  y = .data$auc,
                  fill = .data$run_key)) +
    geom_bar(stat = "identity") +
    xlab("networks") +
    ylab("AUROC") +
    coord_flip() +
    theme(legend.position = "none")

  # AUROC Heatmap
  auroc_heat <- auroc_tibble %>% get_auroc_heat()

  # Extract AU PRROC
  prauc_tibble <- .res_tibble %>%
    unnest(.data$prc) %>%
    select(.data$set_name, .data$bench_name, .data$filter_crit,
           .data$statistic, .data$auc) %>%
    distinct()

  # AU PR Heatmap
  pr_heat <- prauc_tibble %>% get_auroc_heat()

  # get computational time info
  comp_time <- .res_tibble %>%
    # get statistic time from activity
    mutate(statistic_time = .data$activity %>%
             map(function(tib)
               tib %>%
                 select(.data$statistic_time) %>%
                 unique)) %>%
    unnest(.data$statistic_time) %>%
    # calculate regulon size
    group_by(.data$set_name, .data$bench_name, .data$filter_crit) %>%
    mutate(regulon_time = sum(.data$statistic_time)) %>%
    ungroup %>%
    select(.data$set_name, .data$bench_name, .data$statistic,
           .data$filter_crit, .data$statistic_time, .data$regulon_time)

  # Join AUROC, PRAUC, Coverage, and Comp time
  roc_cov <- roc %>%
    group_by(.data$name_lvl) %>%
    summarise(source_cov = .data$coverage,
              condition_cov = n,
              roc_neg=.data$cn) %>%
    distinct() %>%
    ungroup() %>%
    separate(col="name_lvl",
             into=c("set_name", "bench_name", "filter_crit"),
             sep="\\.")

  pr_cov <- pr %>%
    group_by(.data$name_lvl) %>%
    summarise(pr_neg=.data$cn) %>%
    distinct() %>%
    ungroup() %>%
    separate(col="name_lvl",
             into=c("set_name", "bench_name", "filter_crit"),
             sep="\\.")


  summary_table <- auroc_tibble %>%
    inner_join(prauc_tibble %>%
                 rename(pr_auc = .data$auc),
               by = c("set_name", "bench_name", "statistic", "filter_crit")) %>%
    inner_join(x=.,
               y=roc_cov,
               by = c("set_name", "bench_name", "filter_crit")) %>%
    inner_join(x=.,
               y=pr_cov,
               by = c("set_name", "bench_name", "filter_crit")) %>%
    inner_join(x=.,
               y=comp_time,
               by = c("set_name", "bench_name", "filter_crit", "statistic")) %>%
    distinct()

  bench_summary <- list(summary_table, roc_plot, pr_plot,
                        auroc_plot, auroc_heat, pr_heat)

  names(bench_summary) <- c("summary_table", "roc_plot", "pr_plot",
                            "auroc_plot", "auroc_heat", "pr_heat")

  return(bench_summary)
}

format_roc <- function(.res_tibble, roc_column){
  .res_tibble %>%
    select(.data$set_name, .data$bench_name, .data$filter_crit,
           .data$statistic, .data[[roc_column]]) %>%
    unnest(.data[[roc_column]]) %>%
    unite("name_lvl", .data$set_name, .data$bench_name,
          .data$filter_crit, remove = FALSE, sep = ".") %>%
    unite("run_key", .data$set_name, .data$bench_name,
          .data$statistic, .data$filter_crit, remove = FALSE)
}

get_auroc_heat <- function(auroc_tibble){
  auroc_tibble %>%
    select(.data$statistic, .data$auc, .data$filter_crit,
           .data$set_name, .data$bench_name) %>%
    unite("name_lvl", .data$set_name, .data$bench_name, .data$filter_crit) %>%
    pivot_wider(names_from = .data$name_lvl, values_from = .data$auc) %>%
    column_to_rownames(var = "statistic")  %>%
    pheatmap(.,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             treeheight_col = 0,
             treeheight_row = 0,
             display_numbers = TRUE,
             silent = TRUE)
}


# Run benchmark
result <- run_benchmark_mod(
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

saveRDS(result,file=results_out)
