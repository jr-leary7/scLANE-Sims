##### GLM functions #####

simulate_scaffold <- function(ref.dataset = NULL,
                              perc.dyn.genes = NULL,
                              n.cells = NULL) {
  # check inputs
  if (is.null(ref.dataset) | is.null(perc.dyn.genes) | is.null(n.cells)) { stop("You're missing vital parameters for simulate_scaffold().") }
  if (perc.dyn.genes <= 0) { stop("% dynamic genes need to be greater than zero.") }
  if (n.cells <= 0) { stop("Number of cells needs to be greater than zero.") }
  # set up simulation parameters
  n_dyn_genes <- ceiling(perc.dyn.genes * nrow(ref.dataset))
  my_knots <- matrix(runif(2 * n_dyn_genes, 0, 1), ncol = 2, nrow = n_dyn_genes)
  my_theta <- matrix(rnorm(5, 5, 5), ncol = 5, nrow = n_dyn_genes)
  dynamic_params <- list(propGenes = perc.dyn.genes,
                         degree = 2,
                         knots = my_knots,
                         theta = my_theta)
  # simulate 10X dataset
  scaffold_params <- estimateScaffoldParameters(sce = ref.dataset,
                                                sceUMI = TRUE,
                                                useUMI = TRUE,
                                                protocol = "droplet",
                                                numCells = n.cells,
                                                popHet = c(1, 1),
                                                useDynamic = dynamic_params)
  sim_data <- simulateScaffold(scaffoldParams = scaffold_params, originalSCE = ref.dataset)
  sim_data <- logNormCounts(sim_data)
  var_decomp <- modelGeneVar(sim_data)
  top2k_hvgs <- getTopHVGs(var_decomp, n = 2000)
  sim_data <- runPCA(sim_data, subset_row = top2k_hvgs)
  reducedDim(sim_data, "PCAsub") <- reducedDim(sim_data, "PCA")[, 1:10, drop = FALSE]
  sim_data <- runUMAP(sim_data, dimred = "PCAsub", n_dimred = 1:10)
  g <- buildSNNGraph(sim_data, use.dimred = "PCAsub", k = 30)
  clusters <- igraph::cluster_louvain(graph = g)$membership
  colLabels(sim_data) <- factor(clusters)
  colData(sim_data) <- colData(sim_data) %>%
                       as.data.frame() %>%
                       dplyr::mutate(cell_time = as.numeric(gsub("Cell_", "", rownames(.))),
                                     cell_time_normed = cell_time / max(cell_time)) %>%
                       DataFrame()
  return(sim_data)
}

run_scLANE <- function(sim.data = NULL,
                       n.iter = 3,
                       param.list = NULL,
                       n.cores = NULL,
                       scLANE.log = TRUE,
                       scLANE.log.iter = 1000) {
  # check inputs
  if (is.null(sim.data) | is.null(param.list) | is.null(n.cores)) { stop("You failed to provide necessary parameters to run_scLANE().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 8)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "testDynamic_results_raw", "testDynamic_results_tidy", "testSlope_results")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  td_res_raw_list <- vector("list", length = n.iter)
  td_res_tidy_list <- vector("list", length = n.iter)
  ts_res_list <- vector("list", length = n.iter)

  # prepare counts matrix & cell-ordering dataframe
  sim_counts <- as.matrix(t(counts(sim.data)))
  sim_counts <- sim_counts[, which(colSums(sim_counts) > 0)]
  pt_df <- colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        gene_stats <- testDynamic(expr.mat = sim_counts,
                                  pt = pt_df,
                                  parallel.exec = TRUE,
                                  n.cores = n.cores,
                                  n.potential.basis.fns = 5,
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
  return(res_list)
}

run_scLANE_reduced <- function(sim.data = NULL,
                               n.genes.sample = 1000,
                               n.iter = 3,
                               param.list = NULL,
                               n.cores = NULL,
                               scLANE.log = TRUE,
                               scLANE.log.iter = 1000) {
  # check inputs
  if (is.null(sim.data) | is.null(param.list) | is.null(n.cores)) { stop("You failed to provide necessary parameters to run_scLANE().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  if (n.genes.sample <= 0) { stop("n.genes.sample HAS to be positive, come on.") }
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 8)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "testDynamic_results_raw", "testDynamic_results_tidy", "testSlope_results")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  td_res_raw_list <- vector("list", length = n.iter)
  td_res_tidy_list <- vector("list", length = n.iter)
  ts_res_list <- vector("list", length = n.iter)

  # prepare subsampled (preserves % dynamic genes) counts matrix & cell-ordering dataframe (no need to set seed b/c targets takes care of reproducibility)
  p_dynamic <- mean(rowData(a)[, 1] == "Dynamic")
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
  sim_counts <- sim_counts[, which(colSums(sim_counts) > 0)]
  pt_df <- colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        gene_stats <- testDynamic(expr.mat = sim_counts,
                                  pt = pt_df,
                                  parallel.exec = TRUE,
                                  n.cores = n.cores,
                                  n.potential.basis.fns = 5,
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
  return(res_list)
}

##### GEE function #####

simulate_scaffold_GEE <- function(ref.dataset = NULL,
                                  perc.dyn.genes = NULL,
                                  n.cells = NULL,
                                  perc.allocation = NULL,
                                  n.subjects = 6) {
  # check inputs
  if (is.null(ref.dataset) | is.null(perc.dyn.genes) | is.null(n.cells)) { stop("You're missing vital parameters for simulate_scaffold().") }
  if (perc.dyn.genes <= 0) { stop("% dynamic genes need to be greater than zero.") }
  if (n.cells <= 0) { stop("Number of cells needs to be greater than zero.") }
  if (is.null(perc.allocation)) { stop("% allocation must be non-NULL.") }
  if (sum(perc.allocation) != 1) { stop("% allocation must add up to one.") }
  if (length(perc.allocation) != n.subjects) { stop("Each subject must have a % sample allocation value.") }
  obj_list <- vector("list", length = n.subjects)
  # set up simulation parameters -- common across subjects
  n_dyn_genes <- ceiling(perc.dyn.genes * nrow(ref.dataset))
  my_knots <- matrix(runif(2 * n_dyn_genes, 0, 1), ncol = 2, nrow = n_dyn_genes)
  my_theta <- matrix(rnorm(5, 5, 5), ncol = 5, nrow = n_dyn_genes)
  dynamic_params <- list(propGenes = perc.dyn.genes,
                         degree = 2,
                         knots = my_knots,
                         theta = my_theta)
  # simulate 10X dataset
  for (s in seq(n.subjects)) {
    subject_n_cells <- n.cells * perc.allocation[s]
    scaffold_params <- estimateScaffoldParameters(sce = ref.dataset,
                                                  sceUMI = TRUE,
                                                  useUMI = TRUE,
                                                  protocol = "droplet",
                                                  numCells = subject_n_cells,
                                                  popHet = c(1, 1),
                                                  useDynamic = dynamic_params)
    obj_list[[s]] <- simulateScaffold(scaffoldParams = scaffold_params, originalSCE = ref.dataset)
  }
  counts_mat <- purrr::map(obj_list, function(x) { counts(x) }) %>%
                purrr::reduce(cbind) %>%
                as.matrix()  # cast to dense matrix (yikes)
  col_names <- c()
  for (i in seq_along(obj_list)) {
    col_names <- c(col_names, paste0("P", i, "_", colnames(obj_list[[i]])))
  }
  row_names <- rownames(obj_list[[1]])  # common to all sim data
  rownames(counts_mat) <- row_names
  colnames(counts_mat) <- col_names
  row_data <- purrr::map(obj_list, function(x) rowData(x)) %>%
              purrr::reduce(cbind) %>%
              as.data.frame()
  colnames(row_data) <- paste0("geneStatus_P", 1:n.subjects)
  row_data <- DataFrame(row_data)
  col_data <- purrr::map(obj_list, function(x) colData(x)) %>%
              purrr::reduce(rbind) %>%
              as.data.frame()
  subj_names <- c()
  cell_time_normed <- c()
  for (i in seq_along(obj_list)) {
    subj_names <- c(subj_names, rep(paste0("P", i), ncol(obj_list[[i]])))
    cell_time_normed <- c(cell_time_normed, 1:ncol(obj_list[[i]]) / ncol(obj_list[[i]]))
  }
  col_data <- dplyr::mutate(col_data,
                            subject = subj_names,
                            cell_time_normed = cell_time_normed)
  rownames(col_data) <- col_names
  col_data <- DataFrame(col_data)
  sim_data <- SingleCellExperiment(list(counts = counts_mat))
  colData(sim_data) <- col_data
  rowData(sim_data) <- row_data
  # process things
  sim_data <- logNormCounts(sim_data)
  var_decomp <- modelGeneVar(sim_data)
  top2k_hvgs <- getTopHVGs(var_decomp, n = 2000)
  sim_data <- runPCA(sim_data, subset_row = top2k_hvgs)
  reducedDim(sim_data, "PCAsub") <- reducedDim(sim_data, "PCA")[, 1:10, drop = FALSE]
  sim_data <- runUMAP(sim_data, dimred = "PCAsub", n_dimred = 1:10)
  g <- buildSNNGraph(sim_data, use.dimred = "PCAsub", k = 30)
  clusters <- igraph::cluster_louvain(graph = g)$membership
  colLabels(sim_data) <- factor(clusters)
  return(sim_data)
}

run_scLANE_GEE <- function(sim.data = NULL,
                           n.iter = 3,
                           param.list = NULL,
                           n.cores = NULL,
                           scLANE.log = TRUE,
                           scLANE.log.iter = 1000) {
  # check inputs
  if (is.null(sim.data) | is.null(param.list) | is.null(n.cores)) { stop("You failed to provide necessary parameters to run_scLANE().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 8)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "testDynamic_results_raw", "testDynamic_results_tidy", "testSlope_results")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  td_res_raw_list <- vector("list", length = n.iter)
  td_res_tidy_list <- vector("list", length = n.iter)
  ts_res_list <- vector("list", length = n.iter)

  # prepare counts matrix & cell-ordering dataframe
  sim_counts <- as.matrix(t(counts(sim.data)))
  sim_counts <- sim_counts[, which(colSums(sim_counts) > 0)]
  pt_df <- colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
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
  return(res_list)
}

run_scLANE_reduced_GEE <- function(sim.data = NULL,
                                   n.genes.sample = 1000,
                                   n.iter = 3,
                                   param.list = NULL,
                                   n.cores = NULL,
                                   scLANE.log = TRUE,
                                   scLANE.log.iter = 1000) {
  # check inputs
  if (is.null(sim.data) | is.null(param.list) | is.null(n.cores)) { stop("You failed to provide necessary parameters to run_scLANE().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 8)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "testDynamic_results_raw", "testDynamic_results_tidy", "testSlope_results")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  td_res_raw_list <- vector("list", length = n.iter)
  td_res_tidy_list <- vector("list", length = n.iter)
  ts_res_list <- vector("list", length = n.iter)
  # prepare subsampled (preserves % dynamic genes) counts matrix & cell-ordering dataframe (no need to set seed b/c targets takes care of reproducibility)
  p_dynamic <- rowData(a) %>%
               as.data.frame() %>% View()
               select(contains("geneStatus")) %>%
               tidyr::pivot_longer(cols = everything(), values_to = "geneStatus") %>%
               summarise(P = mean(geneStatus == "Dynamic")) %>%
               pull(P)
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
  # prepare counts matrix & cell-ordering dataframe
  sim_counts <- as.matrix(t(counts(sim.data)))
  sim_counts <- sim_counts[, which(colSums(sim_counts) > 0)]
  pt_df <- colData(sim.data) %>%
    as.data.frame() %>%
    dplyr::select(cell_time_normed) %>%
    dplyr::rename(PT = cell_time_normed)
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
  return(res_list)
}
