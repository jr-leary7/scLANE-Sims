### NOT READY YET ###

# run run_monocle3 on subsampled single-subject data
run_monocle3 <- function(sim.data = NULL,
                         n.genes.sample = 1000,
                         n.iter = 1,
                         param.list = NULL,
                         n.cores = 4) {
  # check inputs
  if (is.null(sim.data) || is.null(param.list)) { stop("You failed to provide necessary parameters to run_monocle3().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  if (n.cores <= 0) { stop("n.cores HAS to be positive, come on.") }
  if (n.genes.sample <= 0) { stop("n.genes.sample HAS to be positive, come on.") }
  
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 8)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "monocle3_results_raw", "monocle3_results_tidy", "RMSE_estimates")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  mono_res_raw_list <- vector("list", length = n.iter)
  mono_res_tidy_list <- vector("list", length = n.iter)
  rmse_list <- vector("list", length = n.iter)
  
  # prepare subsampled (preserves % dynamic genes) counts matrix & cell-ordering dataframe -- {targets} takes care of seed
  p_dynamic <- mean(SummarizedExperiment::rowData(sim.data)[, 1] == "Dynamic")
  n_dyn_genes <- ceiling(p_dynamic * n.genes.sample)
  n_norm_genes <- n.genes.sample - n_dyn_genes
  samp_dyn_genes <- SummarizedExperiment::rowData(sim.data) %>%
                    as.data.frame() %>%
                    dplyr::filter(geneStatus == "Dynamic") %>%
                    dplyr::slice_sample(n = n_dyn_genes) %>%
                    rownames(.)
  samp_norm_genes <- SummarizedExperiment::rowData(sim.data) %>%
                     as.data.frame() %>%
                     dplyr::filter(geneStatus == "NotDynamic") %>%
                     dplyr::slice_sample(n = n_norm_genes) %>%
                     rownames(.)
  samp_genes <- c(samp_dyn_genes, samp_norm_genes)
  sim.data <- sim.data[rownames(sim.data) %in% samp_genes, ]
  cds <- monocle3::new_cell_data_set(expression_data = sim.data@assays@data$counts, 
                                     cell_metadata = as.data.frame(SummarizedExperiment::colData(sim.data)), 
                                     gene_metadata = as.data.frame(SummarizedExperiment::rowData(sim.data)))
  
  # run monocle3
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        gene_stats <- monocle3::fit_models(cds,
                                           model_formula_str = "~cell_time_normed",
                                           expression_family = "negbinomial", 
                                           cores = n.cores)
        monocle_coefs <- monocle3::coefficient_table(gene_stats) %>% 
                         dplyr::filter(term == "cell_time_normed") %>% 
                         dplyr::arrange(p_value) %>% 
                         dplyr::mutate(pvalue_adj = stats::p.adjust(p_value, method = "holm"), 
                                       gene_dynamic_overall = dplyr::if_else(pvalue_adj < 0.01, 1, 0))
        end_time <- Sys.time()
        time_diff <- end_time - start_time
      }
    )
    start_time_list[[i]] <- start_time
    end_time_list[[i]] <- end_time
    time_diff_list[[i]] <- time_diff
    mono_res_raw_list[[i]] <- gene_stats
    mono_res_tidy_list[[i]] <- monocle_coefs
    param_list[[i]] <- param.list
    
    ### still need to write helper function to get fitted values from monocle3 results
    
    # rmse_list[[i]] <- purrr::map2(rownames(sim_counts), 
    #                               gene_stats, 
    #                               function(x, y) {
    #                                 rmse <- try({ 
    #                                   yardstick::rmse_vec(truth = sim_counts[x, ], estimate = y$fitted.values)
    #                                 }, silent = TRUE) 
    #                               })
    # names(rmse_list[[i]]) <- rownames(sim_counts)
  }
  
  # set up results list & return
  res_list$sim_parameters <- param_list
  res_list$start_time <- start_time_list
  res_list$end_time <- end_time_list
  res_list$time_diff <- time_diff_list
  res_list$mem_usage <- mem_usage_list
  res_list$monocle3_results_raw <- mono_res_raw_list
  res_list$monocle3_results_tidy <- mono_res_tidy_list
  res_list$RMSE_estimates <- rmse_list
  return(res_list)
}
