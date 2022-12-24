##### Single subject functions #####


### Single-subject data simulation 


simulate_scaffold <- function(ref.dataset = NULL,
                              perc.dyn.genes = NULL,
                              n.cells = NULL) {
  # check inputs
  if (is.null(ref.dataset) | is.null(perc.dyn.genes) | is.null(n.cells)) { stop("You're missing vital parameters for simulate_scaffold().") }
  if (perc.dyn.genes <= 0) { stop("% dynamic genes need to be greater than zero.") }
  if (n.cells <= 0) { stop("Number of cells needs to be greater than zero.") }
  
  # set up simulation parameters
  n_dyn_genes <- ceiling(perc.dyn.genes * nrow(ref.dataset))
  n_possible_dyn_genes <- ceiling((perc.dyn.genes / 0.8) * nrow(ref.dataset))  # make dynamic genes an 80% sample of the total pool of possible dynamic genes
  Q75 <- unname(quantile(rowMeans(counts(ref.dataset)), 0.75))
  high_exp_genes <- rownames(ref.dataset)[rowMeans(SingleCellExperiment::counts(ref.dataset)) > Q75]
  if (length(high_exp_genes) < n_possible_dyn_genes) {
    Q60 <- unname(quantile(rowMeans(SingleCellExperiment::counts(ref.dataset)), 0.6))
    high_exp_genes <- rownames(ref.dataset)[rowMeans(SingleCellExperiment::counts(ref.dataset)) > Q60]
    if (length(high_exp_genes) < n_possible_dyn_genes) {
      stop("Your dataset has too few highly expressed genes to support the number of dynamic genes you want. Please reduce the % dynamic genes parameter.")
    }
  }
  possible_dyn_genes <- sample(high_exp_genes,  # make sure pool of possible dynamic genes has high expression
                               size = n_possible_dyn_genes, 
                               replace = FALSE)
  dyn_genes <- sample(possible_dyn_genes, n_dyn_genes, replace = FALSE)
  my_knots <- matrix(runif(2 * n_dyn_genes, 0, 1), ncol = 2, nrow = n_dyn_genes)
  my_theta <- matrix(rnorm(5, 5, 5), ncol = 5, nrow = n_dyn_genes)
  dynamic_params <- list(propGenes = perc.dyn.genes,
                         dynGenes = dyn_genes, 
                         degree = 2,
                         knots = my_knots,
                         theta = my_theta)
  
  # simulate 10X dataset
  scaffold_params <- scaffold::estimateScaffoldParameters(sce = ref.dataset,
                                                          sceUMI = TRUE,
                                                          useUMI = TRUE,
                                                          protocol = "droplet",
                                                          numCells = n.cells,
                                                          popHet = c(1, 1),
                                                          useDynamic = dynamic_params)
  sim_data <- scaffold::simulateScaffold(scaffoldParams = scaffold_params, originalSCE = ref.dataset)
  
  # typical scran + scater pre-processing pipeline
  sim_data <- sim_data[rowSums(SingleCellExperiment::counts(sim_data) > 0) > 0, ]  # only non-zero genes -- shouldn't drop many dynamic genes thanks to biased selection
  sim_data <- scater::logNormCounts(sim_data)
  var_decomp <- scran::modelGeneVar(sim_data)
  top2k_hvgs <- scran::getTopHVGs(var_decomp, n = 2000)
  sim_data <- scater::runPCA(sim_data, subset_row = top2k_hvgs)
  SingleCellExperiment::reducedDim(sim_data, "PCAsub") <- SingleCellExperiment::reducedDim(sim_data, "PCA")[, 1:10, drop = FALSE]
  sim_data <- scater::runUMAP(sim_data, dimred = "PCAsub", n_dimred = 1:10)
  g <- scran::buildSNNGraph(sim_data, use.dimred = "PCAsub", k = 30)
  clusters <- igraph::cluster_louvain(graph = g)$membership
  SingleCellExperiment::colLabels(sim_data) <- factor(clusters)
  colData(sim_data) <- colData(sim_data) %>%
                       as.data.frame() %>%
                       dplyr::mutate(cell_time = as.numeric(gsub("Cell_", "", rownames(.))),
                                     cell_time_normed = cell_time / max(cell_time)) %>%
                       S4Vectors::DataFrame()
  return(sim_data)
}


### Single-subject scLANE GLM functions


run_scLANE <- function(sim.data = NULL,
                       n.iter = 3,
                       param.list = NULL,
                       n.cores = NULL,
                       scLANE.log = FALSE,
                       scLANE.log.iter = 1000) {
  # check inputs
  if (is.null(sim.data) | is.null(param.list) | is.null(n.cores)) { stop("You failed to provide necessary parameters to run_scLANE().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 9)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "testDynamic_results_raw", "testDynamic_results_tidy", "testSlope_results", "RMSE_estimates")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  td_res_raw_list <- vector("list", length = n.iter)
  td_res_tidy_list <- vector("list", length = n.iter)
  ts_res_list <- vector("list", length = n.iter)
  rmse_list <- vector("list", length = n.iter)

  # prepare counts matrix & cell-ordering dataframe
  sim_counts <- as.matrix(t(counts(sim.data)))
  pt_df <- colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
  
  # run scLANE
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        gene_stats <- testDynamic(expr.mat = sim_counts,
                                  pt = pt_df,
                                  parallel.exec = TRUE,
                                  n.cores = n.cores,
                                  n.potential.basis.fns = 5,
                                  approx.knot = TRUE, 
                                  track.time = TRUE,
                                  log.file = scLANE.log,
                                  log.iter = scLANE.log.iter)
        global_test_results <- getResultsDE(gene_stats, p.adj.method = "bonferroni", fdr.cutoff = 0.01) %>%
                               dplyr::inner_join((rowData(sim.data) %>%
                               as.data.frame() %>%
                               dplyr::mutate(gene = rownames(.))), by = c("Gene" = "gene"))
        slope_test_results <- testSlope(test.dyn.results = gene_stats,
                                        p.adj.method = "bonferroni",
                                        fdr.cutoff = 0.01) %>%
                              dplyr::inner_join((rowData(sim.data) %>%
                              as.data.frame() %>%
                              dplyr::mutate(gene = rownames(.))), by = c("Gene" = "gene"))
        end_time <- Sys.time()
        time_diff <- end_time - start_time
      }
    )
    start_time_list[[i]] <- start_time
    end_time_list[[i]] <- end_time
    time_diff_list[[i]] <- time_diff
    td_res_raw_list[[i]] <- gene_stats
    td_res_tidy_list[[i]] <- global_test_results
    ts_res_list[[i]] <- slope_test_results
    param_list[[i]] <- param.list
    rmse_list[[i]] <- purrr::map2(colnames(sim_counts), 
                                  gene_stats, 
                                  function(x, y) {
                                    rmse <- try({ 
                                      yardstick::rmse_vec(truth = sim_counts[, x], estimate = exp(y$Lineage_A$MARGE_Preds$marge_link_fit))
                                    }, silent = TRUE) 
                                  })
    names(rmse_list[[i]]) <- colnames(sim_counts)
  }
  
  # set up results list & return
  res_list$sim_parameters <- param_list
  res_list$start_time <- start_time_list
  res_list$end_time <- end_time_list
  res_list$time_diff <- time_diff_list
  res_list$mem_usage <- mem_usage_list
  res_list$testDynamic_results_raw <- td_res_raw_list
  res_list$testDynamic_results_tidy <- td_res_tidy_list
  res_list$testSlope_results <- ts_res_list
  res_list$RMSE_estimates <- rmse_list
  return(res_list)
}

run_scLANE_reduced <- function(sim.data = NULL,
                               n.genes.sample = 1000,
                               n.iter = 3,
                               param.list = NULL,
                               n.cores = NULL,
                               scLANE.log = FALSE,
                               scLANE.log.iter = 1000) {
  # check inputs
  if (is.null(sim.data) | is.null(param.list) | is.null(n.cores)) { stop("You failed to provide necessary parameters to run_scLANE().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  if (n.genes.sample <= 0) { stop("n.genes.sample HAS to be positive, come on.") }
  
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 9)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "testDynamic_results_raw", "testDynamic_results_tidy", "testSlope_results", "RMSE_estimates")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  td_res_raw_list <- vector("list", length = n.iter)
  td_res_tidy_list <- vector("list", length = n.iter)
  ts_res_list <- vector("list", length = n.iter)
  rmse_list <- vector("list", length = n.iter)

  # prepare subsampled (preserves % dynamic genes) counts matrix & cell-ordering dataframe (no need to set seed b/c targets takes care of reproducibility)
  p_dynamic <- mean(rowData(sim.data)[, 1] == "Dynamic")
  n_dyn_genes <- ceiling(p_dynamic * n.genes.sample)
  n_norm_genes <- n.genes.sample - n_dyn_genes
  samp_dyn_genes <- rowData(sim.data) %>%
                    as.data.frame() %>%
                    dplyr::filter(geneStatus == "Dynamic") %>%
                    dplyr::slice_sample(n = n_dyn_genes) %>%
                    rownames(.)
  samp_norm_genes <- rowData(sim.data) %>%
                     as.data.frame() %>%
                     dplyr::filter(geneStatus == "NotDynamic") %>%
                     dplyr::slice_sample(n = n_norm_genes) %>%
                     rownames(.)
  samp_genes <- c(samp_dyn_genes, samp_norm_genes)
  sim.data <- sim.data[rownames(sim.data) %in% samp_genes, ]
  sim_counts <- as.matrix(t(SingleCellExperiment::counts(sim.data)))
  pt_df <- colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
  
  # run scLANE
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        gene_stats <- testDynamic(expr.mat = sim_counts,
                                  pt = pt_df,
                                  parallel.exec = TRUE,
                                  n.cores = n.cores,
                                  n.potential.basis.fns = 5,
                                  approx.knot = TRUE, 
                                  track.time = TRUE,
                                  log.file = scLANE.log,
                                  log.iter = scLANE.log.iter)
        global_test_results <- getResultsDE(gene_stats, p.adj.method = "bonferroni", fdr.cutoff = 0.01) %>%
          dplyr::inner_join((rowData(sim.data) %>%
                               as.data.frame() %>%
                               dplyr::mutate(gene = rownames(.))), by = c("Gene" = "gene"))
        slope_test_results <- testSlope(test.dyn.results = gene_stats,
                                        p.adj.method = "bonferroni",
                                        fdr.cutoff = 0.01) %>%
          dplyr::inner_join((rowData(sim.data) %>%
                             as.data.frame() %>%
                             dplyr::mutate(gene = rownames(.))), by = c("Gene" = "gene"))
        end_time <- Sys.time()
        time_diff <- end_time - start_time
      }
    )
    start_time_list[[i]] <- start_time
    end_time_list[[i]] <- end_time
    time_diff_list[[i]] <- time_diff
    td_res_raw_list[[i]] <- gene_stats
    td_res_tidy_list[[i]] <- global_test_results
    ts_res_list[[i]] <- slope_test_results
    param_list[[i]] <- param.list
    rmse_list[[i]] <- purrr::map2(colnames(sim_counts), 
                                  gene_stats, 
                                  function(x, y) {
                                    rmse <- try({ 
                                      yardstick::rmse_vec(truth = sim_counts[, x], estimate = exp(y$Lineage_A$MARGE_Preds$marge_link_fit))
                                    }, silent = TRUE) 
                                  })
    names(rmse_list[[i]]) <- colnames(sim_counts)
  }
  
  # set up results list & return
  res_list$sim_parameters <- param_list
  res_list$start_time <- start_time_list
  res_list$end_time <- end_time_list
  res_list$time_diff <- time_diff_list
  res_list$mem_usage <- mem_usage_list
  res_list$testDynamic_results_raw <- td_res_raw_list
  res_list$testDynamic_results_tidy <- td_res_tidy_list
  res_list$testSlope_results <- ts_res_list
  res_list$RMSE_estimates <- rmse_list
  return(res_list)
}


### Single-subject tradeSeq function


run_tradeSeq_reduced <- function(sim.data = NULL,
                                 n.genes.sample = 1000,
                                 n.iter = 3,
                                 param.list = NULL,
                                 n.cores = NULL) {
  # check inputs
  if (is.null(sim.data) | is.null(param.list) | is.null(n.cores)) { stop("You failed to provide necessary parameters to run_tradeSeq_reduced().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  if (n.genes.sample <= 0) { stop("n.genes.sample HAS to be positive, come on.") }
  
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 8)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "tradeSeq_results_raw", "tradeSeq_results_tidy", "RMSE_estimates")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  ts_res_raw_list <- vector("list", length = n.iter)
  ts_res_tidy_list <- vector("list", length = n.iter)
  rmse_list <- vector("list", length = n.iter)
  
  # prepare subsampled (preserves % dynamic genes) counts matrix & cell-ordering dataframe (no need to set seed b/c targets takes care of reproducibility)
  p_dynamic <- mean(rowData(sim.data)[, 1] == "Dynamic")
  n_dyn_genes <- ceiling(p_dynamic * n.genes.sample)
  n_norm_genes <- n.genes.sample - n_dyn_genes
  samp_dyn_genes <- rowData(sim.data) %>%
                    as.data.frame() %>%
                    dplyr::filter(geneStatus == "Dynamic") %>%
                    dplyr::slice_sample(n = n_dyn_genes) %>%
                    rownames(.)
  samp_norm_genes <- rowData(sim.data) %>%
                     as.data.frame() %>%
                     dplyr::filter(geneStatus == "NotDynamic") %>%
                     dplyr::slice_sample(n = n_norm_genes) %>%
                     rownames(.)
  samp_genes <- c(samp_dyn_genes, samp_norm_genes)
  sim.data <- sim.data[rownames(sim.data) %in% samp_genes, ]
  sim_counts <- as.matrix(t(counts(sim.data)))
  pt_df <- colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
  
  # run tradeSeq
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        doParallel::registerDoParallel(n.cores)
        BiocParallel::register(BiocParallel::DoparParam())
        gene_stats <- tradeSeq::fitGAM(counts = t(sim_counts), 
                                       pseudotime = pt_df, 
                                       cellWeights = matrix(rep(1, nrow(pt_df)), ncol = 1), 
                                       nknots = 6, 
                                       sce = FALSE, 
                                       parallel = TRUE, 
                                       verbose = FALSE, 
                                       BPPARAM = BiocParallel::DoparParam())
        global_test_results <- tradeSeq::associationTest(models = gene_stats, global = TRUE) %>% 
          dplyr::arrange(pvalue) %>% 
          dplyr::mutate(gene = rownames(.), 
                        pvalue_adj = p.adjust(pvalue, method = "bonferroni"), 
                        gene_dynamic_overall = dplyr::case_when(pvalue_adj < 0.01 ~ 1, TRUE ~ 0)) %>% 
          dplyr::relocate(gene) %>% 
          dplyr::inner_join((rowData(sim.data) %>%
                             as.data.frame() %>%
                             dplyr::mutate(gene = rownames(.))), by = "gene")
        end_time <- Sys.time()
        time_diff <- end_time - start_time
      }
    )
    start_time_list[[i]] <- start_time
    end_time_list[[i]] <- end_time
    time_diff_list[[i]] <- time_diff
    ts_res_raw_list[[i]] <- gene_stats
    ts_res_tidy_list[[i]] <- global_test_results
    param_list[[i]] <- param.list
    rmse_list[[i]] <- purrr::map2(colnames(sim_counts), 
                                  gene_stats, 
                                  function(x, y) {
                                    rmse <- try({ 
                                      yardstick::rmse_vec(truth = sim_counts[, x], estimate = y$fitted.values)
                                    }, silent = TRUE) 
                                  })
    names(rmse_list[[i]]) <- colnames(sim_counts)
  }
  
  # set up results list & return
  res_list$sim_parameters <- param_list
  res_list$start_time <- start_time_list
  res_list$end_time <- end_time_list
  res_list$time_diff <- time_diff_list
  res_list$mem_usage <- mem_usage_list
  res_list$tradeSeq_results_raw <- ts_res_raw_list
  res_list$tradeSeq_results_tidy <- ts_res_tidy_list
  res_list$RMSE_estimates <- rmse_list
  return(res_list)
}


##### Multi-subject functions #####


### Multi-subject data simulation


simulate_scaffold_GEE <- function(ref.dataset = NULL,
                                  perc.dyn.genes = NULL,
                                  n.cells = NULL,
                                  perc.allocation = NULL,
                                  n.subjects = 6, 
                                  gene.dyn.threshold = 4) {
  # check inputs
  if (is.null(ref.dataset) | is.null(perc.dyn.genes) | is.null(n.cells)) { stop("You're missing vital parameters for simulate_scaffold().") }
  if (perc.dyn.genes <= 0) { stop("% dynamic genes need to be greater than zero.") }
  if (n.cells <= 0) { stop("Number of cells needs to be greater than zero.") }
  if (is.null(perc.allocation)) { stop("% allocation must be non-NULL.") }
  if (length(perc.allocation) != n.subjects) { stop("Each subject must have a % sample allocation value.") }
  
  # set up simulation parameters -- common across subjects
  obj_list <- vector("list", length = n.subjects)
  n_dyn_genes <- ceiling(perc.dyn.genes * nrow(ref.dataset))
  n_possible_dyn_genes <- ceiling((perc.dyn.genes / 0.8) * nrow(ref.dataset))  # make dynamic genes an 80% sample of the total pool of possible dynamic genes
  Q75 <- unname(quantile(rowMeans(SingleCellExperiment::counts(ref.dataset)), 0.75))
  high_exp_genes <- rownames(ref.dataset)[rowMeans(SingleCellExperiment::counts(ref.dataset)) > Q75]
  if (length(high_exp_genes) < n_possible_dyn_genes) {
    Q60 <- unname(quantile(rowMeans(SingleCellExperiment::counts(ref.dataset)), 0.6))
    high_exp_genes <- rownames(ref.dataset)[rowMeans(SingleCellExperiment::counts(ref.dataset)) > Q60]
    if (length(high_exp_genes) < n_possible_dyn_genes) {
      stop("Your dataset has too few highly expressed genes to support the number of dynamic genes you want. Please reduce the % dynamic genes parameter.")
    }
  }
  # make sure pool of possible dynamic genes has high expression
  possible_dyn_genes <- sample(high_exp_genes,  
                               size = n_possible_dyn_genes, 
                               replace = FALSE)
  # simulate a given dynamic for each possible dynamic gene, which will be preserved across subject IFF the gene is chosen to be dynamic in each subject
  my_knots <- matrix(runif(2 * n_possible_dyn_genes, 0, 1), 
                     ncol = 2, 
                     nrow = n_possible_dyn_genes)
  my_knots <- as.data.frame(my_knots) %>% 
              dplyr::mutate(gene = possible_dyn_genes)
  my_theta <- matrix(rnorm(5, 5, 5), 
                     ncol = 5, 
                     nrow = n_possible_dyn_genes)
  my_theta <- as.data.frame(my_theta) %>% 
              dplyr::mutate(gene = possible_dyn_genes)
  
  # simulate 10X datasets
  for (s in seq(n.subjects)) {
    subject_n_cells <- ceiling(n.cells * perc.allocation[s])
    dyn_genes <- sample(possible_dyn_genes, 
                        n_dyn_genes, 
                        replace = FALSE)
    # fetch knots & theta values for each chosen dynamic gene from the overall pool of possible dynamic gene trends
    dyn_gene_knots <- my_knots %>% 
                      dplyr::filter(gene %in% dyn_genes) %>% 
                      dplyr::select(-gene) %>% 
                      as.matrix()
    dyn_gene_theta <- my_theta %>% 
                      dplyr::filter(gene %in% dyn_genes) %>% 
                      dplyr::select(-gene) %>% 
                      as.matrix()
    dynamic_params <- list(propGenes = perc.dyn.genes,
                           dynGenes = dyn_genes, 
                           degree = 2,
                           knots = dyn_gene_knots,
                           theta = dyn_gene_theta)
    scaffold_params <- estimateScaffoldParameters(sce = ref.dataset,
                                                  sceUMI = TRUE,
                                                  useUMI = TRUE,
                                                  protocol = "droplet",
                                                  numCells = subject_n_cells,
                                                  numGenes = nrow(ref.dataset), 
                                                  popHet = c(1, 1),
                                                  useDynamic = dynamic_params)
    sim_data <- simulateScaffold(scaffoldParams = scaffold_params, originalSCE = ref.dataset)
    obj_list[[s]] <- sim_data
  }
  
  # combine & clean datasets / metadata
  counts_mat <- purrr::map(obj_list, \(x) { SingleCellExperiment::counts(x) }) %>%
                purrr::reduce(cbind) %>%
                as.matrix()  # cast to dense matrix (yikes)
  col_names <- c()
  for (i in seq_along(obj_list)) {
    col_names <- c(col_names, paste0("P", i, "_", colnames(obj_list[[i]])))
  }
  row_names <- rownames(obj_list[[1]])  # common to all sim data
  rownames(counts_mat) <- row_names
  colnames(counts_mat) <- col_names
  row_data <- purrr::map(obj_list, \(x) rowData(x)) %>%
              purrr::reduce(cbind) %>%
              as.data.frame()
  colnames(row_data) <- paste0("geneStatus_P", 1:n.subjects)
  row_data <- dplyr::mutate(row_data, 
                            geneDynamic_n = rowSums(dplyr::across(dplyr::contains("geneStatus_"), \(x) x == "Dynamic")), 
                            geneStatus_overall = dplyr::if_else(geneDynamic_n >= gene.dyn.threshold, "Dynamic", "NotDynamic")) 
  row_data <- S4Vectors::DataFrame(row_data)
  col_data <- purrr::map(obj_list, \(x) colData(x)) %>%
              purrr::reduce(rbind) %>%
              as.data.frame()
  subj_names <- c()
  cell_time_normed <- c()
  for (i in seq_along(obj_list)) {
    subj_names <- c(subj_names, rep(paste0("P", i), ncol(obj_list[[i]])))
    cell_time_normed <- c(cell_time_normed, 1:ncol(obj_list[[i]]) / ncol(obj_list[[i]]))  # cells are ordered by creation time within each object
  }
  col_data <- dplyr::mutate(col_data,
                            subject = subj_names,
                            cell_time_normed = cell_time_normed)
  rownames(col_data) <- col_names
  col_data <- S4Vectors::DataFrame(col_data)
  sim_data <- SingleCellExperiment::SingleCellExperiment(list(counts = counts_mat))
  colData(sim_data) <- col_data
  rowData(sim_data) <- row_data
  
  # process data w/ typical pipeline
  sim_data <- sim_data[rowSums(SingleCellExperiment::counts(sim_data) > 0) > 0, ]  # only non-zero genes -- shouldn't drop many dynamic genes thanks to biased selection
  sim_data <- scater::logNormCounts(sim_data)
  var_decomp <- scran::modelGeneVar(sim_data)
  top2k_hvgs <- scran::getTopHVGs(var_decomp, n = 2000)
  sim_data <- scater::runPCA(sim_data, subset_row = top2k_hvgs)
  SingleCellExperiment::reducedDim(sim_data, "PCAsub") <- SingleCellExperiment::reducedDim(sim_data, "PCA")[, 1:10, drop = FALSE]
  sim_data <- scater::runUMAP(sim_data, dimred = "PCAsub", n_dimred = 1:10)
  g <- scran::buildSNNGraph(sim_data, use.dimred = "PCAsub", k = 30)
  clusters <- igraph::cluster_louvain(graph = g)$membership
  SingleCellExperiment::colLabels(sim_data) <- factor(clusters)
  return(sim_data)
}


### scLANE GEE functions


run_scLANE_GEE <- function(sim.data = NULL,
                           n.iter = 3,
                           param.list = NULL,
                           n.cores = NULL,
                           scLANE.log = FALSE,
                           scLANE.log.iter = 1000) {
  # check inputs
  if (is.null(sim.data) | is.null(param.list) | is.null(n.cores)) { stop("You failed to provide necessary parameters to run_scLANE().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 9)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "testDynamic_results_raw", "testDynamic_results_tidy", "testSlope_results", "RMSE_estimates")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  td_res_raw_list <- vector("list", length = n.iter)
  td_res_tidy_list <- vector("list", length = n.iter)
  ts_res_list <- vector("list", length = n.iter)
  rmse_list <- vector("list", length = n.iter)

  # prepare counts matrix & cell-ordering dataframe
  sim_counts <- as.matrix(t(SingleCellExperiment::counts(sim.data)))
  pt_df <- colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
  
  # run scLANE
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        gene_stats <- testDynamic(expr.mat = sim_counts,
                                  pt = pt_df,
                                  parallel.exec = TRUE,
                                  n.cores = n.cores,
                                  n.potential.basis.fns = 5,
                                  is.gee = TRUE,
                                  id.vec = sim.data$subject,
                                  cor.structure = "exchangeable",
                                  approx.knot = TRUE, 
                                  track.time = TRUE,
                                  log.file = scLANE.log,
                                  log.iter = scLANE.log.iter)
        global_test_results <- getResultsDE(gene_stats, p.adj.method = "bonferroni", fdr.cutoff = 0.01) %>%
                               dplyr::inner_join((rowData(sim.data) %>%
                                                  as.data.frame() %>%
                                                  dplyr::select(geneStatus_overall) %>% 
                                                  dplyr::mutate(gene = rownames(.))), by = c("Gene" = "gene"))
        slope_test_results <- testSlope(test.dyn.results = gene_stats,
                                        p.adj.method = "bonferroni",
                                        fdr.cutoff = 0.01) %>%
                              dplyr::inner_join((rowData(sim.data) %>%
                                                 as.data.frame() %>%
                                                 dplyr::select(geneStatus_overall) %>% 
                                                 dplyr::mutate(gene = rownames(.))), by = c("Gene" = "gene"))
        end_time <- Sys.time()
        time_diff <- end_time - start_time
      }
    )
    start_time_list[[i]] <- start_time
    end_time_list[[i]] <- end_time
    time_diff_list[[i]] <- time_diff
    td_res_raw_list[[i]] <- gene_stats
    td_res_tidy_list[[i]] <- global_test_results
    ts_res_list[[i]] <- slope_test_results
    param_list[[i]] <- param.list
    rmse_list[[i]] <- purrr::map2(colnames(sim_counts), 
                                  gene_stats, 
                                  function(x, y) {
                                    rmse <- try({ 
                                      yardstick::rmse_vec(truth = sim_counts[, x], estimate = exp(y$Lineage_A$MARGE_Preds$marge_link_fit))
                                    }, silent = TRUE) 
                                  })
    names(rmse_list[[i]]) <- colnames(sim_counts)
  }
  
  # set up results list & return
  res_list$sim_parameters <- param_list
  res_list$start_time <- start_time_list
  res_list$end_time <- end_time_list
  res_list$time_diff <- time_diff_list
  res_list$mem_usage <- mem_usage_list
  res_list$testDynamic_results_raw <- td_res_raw_list
  res_list$testDynamic_results_tidy <- td_res_tidy_list
  res_list$testSlope_results <- ts_res_list
  res_list$RMSE_estimates <- rmse_list
  return(res_list)
}

run_scLANE_reduced_GEE <- function(sim.data = NULL,
                                   n.genes.sample = 1000,
                                   n.iter = 3,
                                   param.list = NULL,
                                   n.cores = NULL,
                                   scLANE.log = FALSE,
                                   scLANE.log.iter = 1000) {
  # check inputs
  if (is.null(sim.data) | is.null(param.list) | is.null(n.cores)) { stop("You failed to provide necessary parameters to run_scLANE().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 9)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "testDynamic_results_raw", "testDynamic_results_tidy", "testSlope_results", "RMSE_estimates")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  td_res_raw_list <- vector("list", length = n.iter)
  td_res_tidy_list <- vector("list", length = n.iter)
  ts_res_list <- vector("list", length = n.iter)
  rmse_list <- vector("list", length = n.iter)
  
  # prepare subsampled (preserves % dynamic genes) counts matrix & cell-ordering dataframe (no need to set seed b/c targets takes care of reproducibility)
  p_dynamic <- rowData(sim.data) %>%
               as.data.frame() %>% 
               dplyr::select(dplyr::contains("geneStatus_P")) %>%
               tidyr::pivot_longer(cols = tidyselect::everything(), values_to = "geneStatus") %>%
               dplyr::summarise(P = mean(geneStatus == "Dynamic")) %>%
               dplyr::pull(P)
  n_dyn_genes <- ceiling(p_dynamic * n.genes.sample)
  n_norm_genes <- n.genes.sample - n_dyn_genes
  samp_dyn_genes <- rowData(sim.data) %>%
                    as.data.frame() %>%
                    dplyr::filter(geneStatus_overall == "Dynamic") %>%
                    dplyr::slice_sample(n = n_dyn_genes) %>%
                    rownames(.)
  samp_norm_genes <- rowData(sim.data) %>%
                     as.data.frame() %>%
                     dplyr::filter(geneStatus_overall == "NotDynamic") %>%
                     dplyr::slice_sample(n = n_norm_genes) %>%
                     rownames(.)
  samp_genes <- c(samp_dyn_genes, samp_norm_genes)
  sim.data <- sim.data[rownames(sim.data) %in% samp_genes, ]
  
  # prepare counts matrix & cell-ordering dataframe
  sim_counts <- as.matrix(t(SingleCellExperiment::counts(sim.data)))
  pt_df <- colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
  
  # run scLANE
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        gene_stats <- testDynamic(expr.mat = sim_counts,
                                  pt = pt_df,
                                  parallel.exec = FALSE,
                                  n.cores = n.cores,
                                  n.potential.basis.fns = 5,
                                  is.gee = TRUE,
                                  id.vec = sim.data$subject,
                                  cor.structure = "exchangeable",
                                  approx.knot = TRUE, 
                                  track.time = TRUE,
                                  log.file = scLANE.log,
                                  log.iter = scLANE.log.iter)
        global_test_results <- getResultsDE(gene_stats, p.adj.method = "bonferroni", fdr.cutoff = 0.01) %>%
                               dplyr::inner_join((rowData(sim.data) %>%
                                                  as.data.frame() %>%
                                                  dplyr::select(geneStatus_overall) %>% 
                                                  dplyr::mutate(gene = rownames(.))), 
                                                 by = c("Gene" = "gene"))
        slope_test_results <- testSlope(test.dyn.results = gene_stats,
                                        p.adj.method = "bonferroni",
                                        fdr.cutoff = 0.01) %>%
                              dplyr::inner_join((rowData(sim.data) %>%
                                                 as.data.frame() %>%
                                                 dplyr::select(geneStatus_overall) %>% 
                                                 dplyr::mutate(gene = rownames(.))), 
                                                by = c("Gene" = "gene"))
        end_time <- Sys.time()
        time_diff <- end_time - start_time
      }
    )
    start_time_list[[i]] <- start_time
    end_time_list[[i]] <- end_time
    time_diff_list[[i]] <- time_diff
    td_res_raw_list[[i]] <- gene_stats
    td_res_tidy_list[[i]] <- global_test_results
    ts_res_list[[i]] <- slope_test_results
    param_list[[i]] <- param.list
    rmse_list[[i]] <- purrr::map2(colnames(sim_counts), 
                                  gene_stats, 
                                  function(x, y) {
                                    rmse <- try({ 
                                      yardstick::rmse_vec(truth = sim_counts[, x], estimate = exp(y$Lineage_A$MARGE_Preds$marge_link_fit))
                                    }, silent = TRUE) 
                                  })
    names(rmse_list[[i]]) <- colnames(sim_counts)
  }
  
  # set up results list & return
  res_list$sim_parameters <- param_list
  res_list$start_time <- start_time_list
  res_list$end_time <- end_time_list
  res_list$time_diff <- time_diff_list
  res_list$mem_usage <- mem_usage_list
  res_list$testDynamic_results_raw <- td_res_raw_list
  res_list$testDynamic_results_tidy <- td_res_tidy_list
  res_list$testSlope_results <- ts_res_list
  res_list$RMSE_estimates <- rmse_list
  return(res_list)
}


### Multi-sample tradeSeq function


run_tradeSeq_reduced_GEE <- function(sim.data = NULL,
                                     n.genes.sample = 1000,
                                     n.iter = 3,
                                     param.list = NULL,
                                     n.cores = NULL) {
  # check inputs
  if (is.null(sim.data) | is.null(param.list) | is.null(n.cores)) { stop("You failed to provide necessary parameters to run_tradeSeq_reduced().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  if (n.genes.sample <= 0) { stop("n.genes.sample HAS to be positive, come on.") }
  
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 8)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "tradeSeq_results_raw", "tradeSeq_results_tidy", "RMSE_estimates")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  ts_res_raw_list <- vector("list", length = n.iter)
  ts_res_tidy_list <- vector("list", length = n.iter)
  rmse_list <- vector("list", length = n.iter)
  
  # prepare subsampled (preserves % dynamic genes) counts matrix & cell-ordering dataframe (no need to set seed b/c targets takes care of reproducibility)
  p_dynamic <- rowData(sim.data) %>%
               as.data.frame() %>% 
               dplyr::select(dplyr::contains("geneStatus_P")) %>%
               tidyr::pivot_longer(cols = tidyselect::everything(), values_to = "geneStatus") %>%
               dplyr::summarise(P = mean(geneStatus == "Dynamic")) %>%
               dplyr::pull(P)
  n_dyn_genes <- ceiling(p_dynamic * n.genes.sample)
  n_norm_genes <- n.genes.sample - n_dyn_genes
  samp_dyn_genes <- rowData(sim.data) %>%
                    as.data.frame() %>%
                    dplyr::filter(geneStatus_overall == "Dynamic") %>%
                    dplyr::slice_sample(n = n_dyn_genes) %>%
                    rownames(.)
  samp_norm_genes <- rowData(sim.data) %>%
                     as.data.frame() %>%
                     dplyr::filter(geneStatus_overall == "NotDynamic") %>%
                     dplyr::slice_sample(n = n_norm_genes) %>%
                     rownames(.)
  samp_genes <- c(samp_dyn_genes, samp_norm_genes)
  sim.data <- sim.data[rownames(sim.data) %in% samp_genes, ]
  
  # prepare counts matrix & cell-ordering dataframe
  sim_counts <- as.matrix(t(SingleCellExperiment::counts(sim.data)))
  pt_df <- colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
  
  # run tradeSeq
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        doParallel::registerDoParallel(n.cores)
        BiocParallel::register(BiocParallel::DoparParam())
        gene_stats <- tradeSeq::fitGAM(counts = t(sim_counts), 
                                       pseudotime = pt_df, 
                                       cellWeights = matrix(rep(1, nrow(pt_df)), ncol = 1), 
                                       nknots = 6, 
                                       sce = FALSE, 
                                       parallel = TRUE, 
                                       verbose = FALSE, 
                                       BPPARAM = BiocParallel::DoparParam())
        global_test_results <- tradeSeq::associationTest(models = gene_stats, global = TRUE) %>% 
                               dplyr::arrange(pvalue) %>% 
                               dplyr::mutate(gene = rownames(.), 
                                             pvalue_adj = p.adjust(pvalue, method = "bonferroni"), 
                                             gene_dynamic_overall = dplyr::case_when(pvalue_adj < 0.01 ~ 1, TRUE ~ 0)) %>% 
                               dplyr::relocate(gene) %>% 
                               dplyr::inner_join((rowData(sim.data) %>%
                                                  as.data.frame() %>%
                                                  dplyr::select(geneStatus_overall) %>% 
                                                  dplyr::mutate(gene = rownames(.))), by = "gene")
        end_time <- Sys.time()
        time_diff <- end_time - start_time
      }
    )
    start_time_list[[i]] <- start_time
    end_time_list[[i]] <- end_time
    time_diff_list[[i]] <- time_diff
    ts_res_raw_list[[i]] <- gene_stats
    ts_res_tidy_list[[i]] <- global_test_results
    param_list[[i]] <- param.list
    rmse_list[[i]] <- purrr::map2(colnames(sim_counts), 
                                  gene_stats, 
                                  function(x, y) {
                                    rmse <- try({ 
                                      yardstick::rmse_vec(truth = sim_counts[, x], estimate = y$fitted.values)
                                    }, silent = TRUE) 
                                  })
    names(rmse_list[[i]]) <- colnames(sim_counts)
  }
  
  # set up results list & return
  res_list$sim_parameters <- param_list
  res_list$start_time <- start_time_list
  res_list$end_time <- end_time_list
  res_list$time_diff <- time_diff_list
  res_list$mem_usage <- mem_usage_list
  res_list$tradeSeq_results_raw <- ts_res_raw_list
  res_list$tradeSeq_results_tidy <- ts_res_tidy_list
  res_list$RMSE_estimates <- rmse_list
  return(res_list)
}


### scLANE GLMM functions 


run_scLANE_GLMM <- function(sim.data = NULL,
                            basis.adaptive = TRUE, 
                            n.iter = 3,
                            param.list = NULL,
                            n.cores = NULL,
                            scLANE.log = FALSE,
                            scLANE.log.iter = 1000) {
  # check inputs
  if (is.null(sim.data) | is.null(param.list) | is.null(n.cores)) { stop("You failed to provide necessary parameters to run_scLANE().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 9)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "testDynamic_results_raw", "testDynamic_results_tidy", "testSlope_results", "RMSE_estimates")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  td_res_raw_list <- vector("list", length = n.iter)
  td_res_tidy_list <- vector("list", length = n.iter)
  ts_res_list <- vector("list", length = n.iter)
  rmse_list <- vector("list", length = n.iter)
  
  # prepare counts matrix & cell-ordering dataframe
  sim_counts <- as.matrix(t(SingleCellExperiment::counts(sim.data)))
  pt_df <- colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
  
  # run scLANE
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        gene_stats <- testDynamic(expr.mat = sim_counts,
                                  pt = pt_df,
                                  parallel.exec = TRUE,
                                  n.cores = n.cores,
                                  n.potential.basis.fns = 3,
                                  is.glmm = TRUE,
                                  id.vec = sim.data$subject, 
                                  glmm.adaptive = basis.adaptive, 
                                  approx.knot = TRUE, 
                                  track.time = TRUE,
                                  log.file = scLANE.log,
                                  log.iter = scLANE.log.iter)
        global_test_results <- getResultsDE(gene_stats, p.adj.method = "bonferroni", fdr.cutoff = 0.01) %>%
                               dplyr::inner_join((rowData(sim.data) %>%
                                                  as.data.frame() %>%
                                                  dplyr::select(geneStatus_overall) %>% 
                                                  dplyr::mutate(gene = rownames(.))), 
                                                 by = c("Gene" = "gene"))
        slope_test_results <- testSlope(test.dyn.results = gene_stats,
                                        p.adj.method = "bonferroni",
                                        fdr.cutoff = 0.01) %>%
                              dplyr::inner_join((rowData(sim.data) %>%
                                                 as.data.frame() %>%
                                                 dplyr::select(geneStatus_overall) %>% 
                                                 dplyr::mutate(gene = rownames(.))), 
                                                by = c("Gene" = "gene"))
        end_time <- Sys.time()
        time_diff <- end_time - start_time
      }
    )
    start_time_list[[i]] <- start_time
    end_time_list[[i]] <- end_time
    time_diff_list[[i]] <- time_diff
    td_res_raw_list[[i]] <- gene_stats
    td_res_tidy_list[[i]] <- global_test_results
    ts_res_list[[i]] <- slope_test_results
    param_list[[i]] <- param.list
    rmse_list[[i]] <- purrr::map2(colnames(sim_counts), 
                                  gene_stats, 
                                  function(x, y) {
                                    rmse <- try({ 
                                      yardstick::rmse_vec(truth = sim_counts[, x], estimate = exp(y$Lineage_A$MARGE_Preds$marge_link_fit))
                                    }, silent = TRUE) 
                                  })
    names(rmse_list[[i]]) <- colnames(sim_counts)
  }
  
  # set up results list & return
  res_list$sim_parameters <- param_list
  res_list$start_time <- start_time_list
  res_list$end_time <- end_time_list
  res_list$time_diff <- time_diff_list
  res_list$mem_usage <- mem_usage_list
  res_list$testDynamic_results_raw <- td_res_raw_list
  res_list$testDynamic_results_tidy <- td_res_tidy_list
  res_list$testSlope_results <- ts_res_list
  res_list$RMSE_estimates <- rmse_list
  return(res_list)
}

run_scLANE_reduced_GLMM <- function(sim.data = NULL,
                                    basis.adaptive = TRUE, 
                                    n.genes.sample = 1000,
                                    n.iter = 3,
                                    param.list = NULL,
                                    n.cores = NULL,
                                    scLANE.log = FALSE,
                                    scLANE.log.iter = 1000) {
  # check inputs
  if (is.null(sim.data) | is.null(param.list) | is.null(n.cores)) { stop("You failed to provide necessary parameters to run_scLANE().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 9)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "testDynamic_results_raw", "testDynamic_results_tidy", "testSlope_results", "RMSE_estimates")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  td_res_raw_list <- vector("list", length = n.iter)
  td_res_tidy_list <- vector("list", length = n.iter)
  ts_res_list <- vector("list", length = n.iter)
  rmse_list <- vector("list", length = n.iter)
  
  # prepare subsampled (preserves % dynamic genes) counts matrix & cell-ordering dataframe (no need to set seed b/c targets takes care of reproducibility)
  p_dynamic <- rowData(sim.data) %>%
               as.data.frame() %>% 
               dplyr::select(dplyr::contains("geneStatus_P")) %>%
               tidyr::pivot_longer(cols = tidyselect::everything(), values_to = "geneStatus") %>%
               dplyr::summarise(P = mean(geneStatus == "Dynamic")) %>%
               dplyr::pull(P)
  n_dyn_genes <- ceiling(p_dynamic * n.genes.sample)
  n_norm_genes <- n.genes.sample - n_dyn_genes
  samp_dyn_genes <- rowData(sim.data) %>%
                    as.data.frame() %>%
                    dplyr::filter(geneStatus_overall == "Dynamic") %>%
                    dplyr::slice_sample(n = n_dyn_genes) %>%
                    rownames(.)
  samp_norm_genes <- rowData(sim.data) %>%
                     as.data.frame() %>%
                     dplyr::filter(geneStatus_overall == "NotDynamic") %>%
                     dplyr::slice_sample(n = n_norm_genes) %>%
                     rownames(.)
  samp_genes <- c(samp_dyn_genes, samp_norm_genes)
  sim.data <- sim.data[rownames(sim.data) %in% samp_genes, ]
  
  # prepare counts matrix & cell-ordering dataframe
  sim_counts <- as.matrix(t(SingleCellExperiment::counts(sim.data)))
  pt_df <- colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
  
  # run scLANE
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        gene_stats <- testDynamic(expr.mat = sim_counts,
                                  pt = pt_df,
                                  parallel.exec = TRUE,
                                  n.cores = n.cores,
                                  n.potential.basis.fns = 3,
                                  is.glmm = TRUE,
                                  id.vec = sim.data$subject, 
                                  glmm.adaptive = basis.adaptive, 
                                  approx.knot = TRUE, 
                                  track.time = TRUE,
                                  log.file = scLANE.log,
                                  log.iter = scLANE.log.iter)
        global_test_results <- getResultsDE(gene_stats, 
                                            p.adj.method = "bonferroni", 
                                            fdr.cutoff = 0.01) %>%
                               dplyr::inner_join((rowData(sim.data) %>%
                                                  as.data.frame() %>%
                                                  dplyr::select(geneStatus_overall) %>% 
                                                  dplyr::mutate(gene = rownames(.))), 
                                                 by = c("Gene" = "gene"))
        slope_test_results <- testSlope(test.dyn.results = gene_stats,
                                        p.adj.method = "bonferroni",
                                        fdr.cutoff = 0.01) %>%
                              dplyr::inner_join((rowData(sim.data) %>%
                                                 as.data.frame() %>%
                                                 dplyr::select(geneStatus_overall) %>% 
                                                 dplyr::mutate(gene = rownames(.))), 
                                                by = c("Gene" = "gene"))
        end_time <- Sys.time()
        time_diff <- end_time - start_time
      }
    )
    start_time_list[[i]] <- start_time
    end_time_list[[i]] <- end_time
    time_diff_list[[i]] <- time_diff
    td_res_raw_list[[i]] <- gene_stats
    td_res_tidy_list[[i]] <- global_test_results
    ts_res_list[[i]] <- slope_test_results
    param_list[[i]] <- param.list
    rmse_list[[i]] <- purrr::map2(colnames(sim_counts), 
                                  gene_stats, 
                                  function(x, y) {
                                    rmse <- try({ 
                                      yardstick::rmse_vec(truth = sim_counts[, x], estimate = exp(y$Lineage_A$MARGE_Preds$marge_link_fit))
                                    }, silent = TRUE) 
                                  })
    names(rmse_list[[i]]) <- colnames(sim_counts)
  }
  
  # set up results list & return
  res_list$sim_parameters <- param_list
  res_list$start_time <- start_time_list
  res_list$end_time <- end_time_list
  res_list$time_diff <- time_diff_list
  res_list$mem_usage <- mem_usage_list
  res_list$testDynamic_results_raw <- td_res_raw_list
  res_list$testDynamic_results_tidy <- td_res_tidy_list
  res_list$testSlope_results <- ts_res_list
  res_list$RMSE_estimates <- rmse_list
  return(res_list)
}


### Lamian functions 


get_pop_fit_Lamian <- function(lamian.res = NULL) {
  # check inputs 
  if (is.null(lamian.res)) { stop("lamian.res must be non-NULL.") }
  
  # get predicted values from GAMs per-gene
  design <- lamian.res$design
  design <- design[, 1, drop = FALSE]
  knotnum <- lamian.res$knotnum
  pt_vec <- lamian.res$pseudotime
  pt_vec <- pt_vec[order(pt_vec)]
  gene <- rownames(lamian.res$statistics)
  fit_list <- purrr::map(gene, function(g) {
    tmp <- matrix(lamian.res$parameter[[g]]$beta, ncol = knotnum[g] + 4)
    beta <- as.vector(tmp[1, ])
    x <- purrr::map(rownames(design), function(i) {
      kronecker(diag(knotnum[g] + 4), design[i, , drop = FALSE])
    })
    if (knotnum[g] == 0) {
      phi <- splines::bs(pt_vec, intercept = TRUE)
    } else {
      knots <- seq(min(pt_vec), max(pt_vec), length.out = knotnum[g] + 2)[2:(knotnum[g] + 1)]
      phi <- splines::bs(pt_vec, knots = knots, intercept = TRUE)
    }
    i <- x[[1]]
    if (ncol(phi) == nrow(i)) {
      fit <- phi %*% i %*% beta
    } else {
      fit <- phi %*% t(i) %*% beta
    }
    fit_df <- data.frame(fit_link = as.numeric(fit))
    rownames(fit_df) <- rownames(fit)
    return(fit_df)
  })
  names(fit_list) <- gene
  return(fit_list)
}

get_preds_Lamian <- function(lamian.res = NULL) {
  # check inputs 
  if (is.null(lamian.res)) { stop("lamian.res must be non-NULL.") }
  
  # run 
  expr <- lamian.res$expr
  design <- lamian.res$design
  cellanno <- lamian.res$cellanno
  knotnum <- lamian.res$knotnum
  pt_vec <- lamian.res$pseudotime[colnames(expr)]
  gene <- rownames(expr)
  
  phi_list <- purrr::map(sort(unique(knotnum)), function(k) {
    if (k == 0) {
      phi <- splines::bs(pt_vec, intercept = TRUE)
    } else {
      knots <- seq(min(pt_vec), max(pt_vec), length.out = k + 2)[2:(k + 1)]
      phi <- splines::bs(pt_vec, knots = knots, intercept = TRUE)
    }
  })
  names(phi_list) <- as.character(sort(unique(knotnum)))
  
  sname <- purrr::map(rownames(design), function(i) {
    cellanno[cellanno[, 2] == i, 1]
  })
  names(sname) <- rownames(design)
  
  k <- unique(knotnum)[2]
  pred <- purrr::map(unique(knotnum), function(k) {
    genesub <- names(knotnum)[knotnum == k]
    B <- purrr::map(genesub, function(g) {
      lamian.res$parameter[[g]]$beta
    }) %>% 
      purrr::reduce(rbind)
    if (length(genesub) == 1) {
      B <- matrix(B, nrow = 1)
    }
    rownames(B) <- genesub
    
    omega <- purrr::map(genesub, function(g) {
      lamian.res$parameter[[g]]$omega
    }) %>% 
      purrr::reduce(rbind)
    if (length(genesub) == 1) {
      omega <- matrix(omega, nrow = 1)
    }
    rownames(omega) <- genesub
    
    phi <- phi_list[[as.character(k)]]
    phi <- purrr::map(rownames(design), function(ss) {
      phi[sname[[ss]], ]
    })
    names(phi) <- rownames(design)
    
    xs <- purrr::map(rownames(design), function(i) {
      kronecker(diag(k + 4), design[i, 1, drop = FALSE])
    })
    names(xs) <- rownames(design)
    
    phi_phi <- purrr::map(rownames(design), function(s) {
      t(phi[[s]]) %*% phi[[s]]
    })
    names(phi_phi) <- rownames(design)
    
    phiX <- purrr::map(rownames(design), function(s) {
      phi[[s]] %*% t(xs[[s]])
    })
    names(phiX) <- rownames(design)
    
    s <- rownames(design)[1]
    predtmp <- purrr::map(rownames(design), function(s) {
      sexpr <- expr[genesub, , drop = FALSE]
      sexpr_phibx <- sexpr[genesub, cellanno[, 2] == s, drop = FALSE] - (B[genesub, ] %*% t(phiX[[s]]))
      nb <- k + 4
      oinv <- purrr::map(genesub, function(g) {
        chol2inv(chol(matrix(omega[g, , drop = FALSE], nrow = nb)))
      })
      names(oinv) <- genesub
      Jchol <- purrr::map(genesub, function(g) {
        chol(phi_phi[[s]] + oinv[[g]])
      })
      names(Jchol) <- genesub
      Jsolve <- purrr::map(genesub, function(g) {
        as.numeric(chol2inv(Jchol[[g]]))
      }) %>% 
        purrr::reduce(cbind)
      if (length(genesub) == 1) {
        Jsolve <- matrix(Jsolve, ncol = 1)
      }
      colnames(Jsolve) <- genesub
      
      K <- tcrossprod(t(phi[[s]]), sexpr_phibx)
      JK <- rowsum((Jsolve * K[rep(seq_len(nb), nb), , drop = FALSE]), rep(seq_len(nb), each =  nb)) 
      res <- t(phi[[s]] %*% JK)
      return(res)
    }) %>% 
      purrr::reduce(cbind)
  }) %>% 
    purrr::reduce(rbind)
  pred <- pred[gene, colnames(expr), drop = FALSE]
  populationFit <- t(get_pop_fit_Lamian(lamian.res = lamian.res) %>% purrr::reduce(cbind))
  pred_res <- pred + populationFit
  return(pred_res)
}

compute_resp_preds_Lamian <- function(lamian.preds, orig.sce) {
  # check inputs 
  if (is.null(lamian.preds) || is.null(orig.sce)) { stop("Inputs must be non-NULL.") }
  size_factors <- BiocGenerics::sizeFactors(orig.sce)
  
  # compute fitted values on response scale (raw counts)
  pred_res_resp <- purrr::map(seq(nrow(lamian.preds)), function(r) {
    fitted_vals <- lamian.preds[r, , drop = FALSE]
    fitted_vals_resp <- t(t(2^(fitted_vals) - 1) * size_factors)
    return(fitted_vals_resp)
  }) %>% 
    purrr::reduce(rbind)
  rownames(pred_res_resp) <- rownames(lamian.preds)
  return(pred_res_resp)
}

run_Lamian_reduced_GEE <- function(sim.data = NULL, 
                                   n.genes.sample = 1000,
                                   n.iter = 3,
                                   param.list = NULL,
                                   n.cores = NULL, 
                                   lamian.iter = 100) {
  # check inputs
  if (is.null(sim.data) | is.null(param.list) | is.null(n.cores)) { stop("You failed to provide necessary parameters to run_tradeSeq_reduced().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  if (n.genes.sample <= 0) { stop("n.genes.sample HAS to be positive, come on.") }
  
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 8)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "Lamian_results_raw", "Lamian_results_tidy", "RMSE_estimates")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  lamian_res_raw_list <- vector("list", length = n.iter)
  lamian_res_tidy_list <- vector("list", length = n.iter)
  rmse_list <- vector("list", length = n.iter)
  
  # prepare subsampled (preserves % dynamic genes) counts matrix & cell-ordering dataframe (no need to set seed b/c targets takes care of reproducibility)
  p_dynamic <- rowData(sim.data) %>%
               as.data.frame() %>% 
               dplyr::select(dplyr::contains("geneStatus_P")) %>%
               tidyr::pivot_longer(cols = tidyselect::everything(), values_to = "geneStatus") %>%
               dplyr::summarise(P = mean(geneStatus == "Dynamic")) %>%
               dplyr::pull(P)
  n_dyn_genes <- ceiling(p_dynamic * n.genes.sample)
  n_norm_genes <- n.genes.sample - n_dyn_genes
  samp_dyn_genes <- rowData(sim.data) %>%
                    as.data.frame() %>%
                    dplyr::filter(geneStatus_overall == "Dynamic") %>%
                    dplyr::slice_sample(n = n_dyn_genes) %>%
                    rownames(.)
  samp_norm_genes <- rowData(sim.data) %>%
                     as.data.frame() %>%
                     dplyr::filter(geneStatus_overall == "NotDynamic") %>%
                     dplyr::slice_sample(n = n_norm_genes) %>%
                     rownames(.)
  samp_genes <- c(samp_dyn_genes, samp_norm_genes)
  sim.data <- sim.data[rownames(sim.data) %in% samp_genes, ]
  
  # prepare counts matrix & cell-ordering dataframe
  cell_anno <- colData(sim.data) %>% 
               as.data.frame() %>% 
               dplyr::mutate(Cell = rownames(.), .before = 1) %>% 
               dplyr::select(Cell, Sample = subject) %>% 
               dplyr::mutate(dplyr::across(tidyselect::everything(), as.character)) %>% 
               magrittr::set_rownames(NULL)
  cell_pt <- as.integer(rank(sim.data$cell_time_normed))
  names(cell_pt) <- cell_anno$Cell
  samp_design <- data.frame(intercept = rep(1, length(unique(cell_anno$Sample)))) %>%
                 magrittr::set_rownames(unique(cell_anno$Sample)) %>%
                 as.matrix()
  sim_counts <- as.matrix(SingleCellExperiment::logcounts(sim.data))
  
  # run Lamian
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        gene_stats <- Lamian::lamian_test(expr = sim_counts,
                                          cellanno = cell_anno,
                                          pseudotime = cell_pt,
                                          design = samp_design,
                                          test.type = "time", 
                                          test.method = "permutation", 
                                          permuiter = lamian.iter, 
                                          ncores = n.cores)
        global_test_results <- gene_stats$statistics %>% 
                               dplyr::rename(pvalue = pval.overall, 
                                             Zstat = z.overall) %>% 
                               dplyr::arrange(pvalue) %>% 
                               dplyr::mutate(pvalue_adj = stats::p.adjust(pvalue, method = "bonferroni")) %>% 
                               dplyr::mutate(gene = rownames(.), 
                                             gene_dynamic_overall = dplyr::case_when(pvalue_adj < 0.01 ~ 1, TRUE ~ 0)) %>% 
                               dplyr::select(gene, 
                                             Zstat, 
                                             pvalue,
                                             pvalue_adj, 
                                             gene_dynamic_overall) %>% 
                               dplyr::inner_join((rowData(sim.data) %>%
                                                  as.data.frame() %>%
                                                  dplyr::select(geneStatus_overall) %>% 
                                                  dplyr::mutate(gene = rownames(.))), 
                                                 by = "gene")
        end_time <- Sys.time()
        time_diff <- end_time - start_time
      }
    )
    gene_stats$fitted_values_link <- get_preds_Lamian(lamian.res = gene_stats)
    gene_stats$fitted_values_resp <- compute_resp_preds_Lamian(lamian.preds = gene_stats$fitted_values_link, orig.sce = sim.data)
    start_time_list[[i]] <- start_time
    end_time_list[[i]] <- end_time
    time_diff_list[[i]] <- time_diff
    lamian_res_raw_list[[i]] <- gene_stats
    lamian_res_tidy_list[[i]] <- global_test_results
    param_list[[i]] <- param.list
    rmse_list[[i]] <- purrr::map(rownames(gene_stats$fitted_values_resp), 
                                 function(x) {
                                   rmse <- try({ 
                                     yardstick::rmse_vec(truth = sim_counts[x, ], estimate = gene_stats$fitted_values_resp[x, ])
                                   }, silent = TRUE) 
                                 })
    names(rmse_list[[i]]) <- rownames(gene_stats$fitted_values_resp)
  }
  
  # set up results list & return
  res_list$sim_parameters <- param_list
  res_list$start_time <- start_time_list
  res_list$end_time <- end_time_list
  res_list$time_diff <- time_diff_list
  res_list$mem_usage <- mem_usage_list
  res_list$Lamian_results_raw <- lamian_res_raw_list
  res_list$Lamian_results_tidy <- lamian_res_tidy_list
  res_list$RMSE_estimates <- rmse_list
  return(res_list)
}
