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
