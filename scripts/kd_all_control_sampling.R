library(decoupleR)
library(decoupleRBench)
library(dplyr)
library(purrr)
library(tibble)
library(yardstick)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(tidyselect)


# This script calculates AUROC and PRC for each method using knocked down samples as
# the test and alll control samples across all experiments as the reference.

args = commandArgs(trailingOnly=TRUE)

bench_results_h = args[1]
downsampling = eval(args[2])
results_out = args[3]

bench_results <- readRDS(bench_results_h)

calc_curve = function(df,
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
    prepare_for_roc(., filter_tn = F)


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


prepare_for_roc = function(df, filter_tn = FALSE) {
  res = df %>%
    dplyr::mutate(response = case_when((.data$source == .data$target & .data$sign ==-1) ~ 1,
                                       (.data$sign == 1) ~ 0 ),
                  predictor = abs(.data$score))
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


res <- bench_results@bench_res %>% mutate(roc = .data$activity %>% map(~calc_curve(df=.x,downsampling=downsampling,times=100,curve="ROC")),
                                      prc = .data$activity %>% map(~calc_curve(df=.x,downsampling=downsampling,times=100,curve="PR")))

bench_result_updated <- methods::new("BenchResult",bench_res=res,summary=res %>% get_bench_summary(),design=bench_results@design)

# Save result
saveRDS(bench_result_updated,file=results_out) 
