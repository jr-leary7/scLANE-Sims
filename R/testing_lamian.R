library(dplyr)
library(Lamian)
library(targets)
library(SingleCellExperiment)

tar_load(lung_sim_DEG_20_CELLS_100_balanced)

tar_load(scLANE_res_lung_DEG_20_CELLS_100_balanced)
scLANE_res_lung_DEG_20_CELLS_100_balanced$testDynamic_results_tidy[[1]] %>% View()

tar_load(tradeSeq_res_lung_DEG_10_CELLS_1000_unbalanced)
tradeSeq_res_lung_DEG_10_CELLS_1000_unbalanced$tradeSeq_results_tidy[[1]] %>% View()


lamian_res <- run_Lamian_reduced_GEE(sim.data = lung_sim_DEG_20_CELLS_100_balanced,
                                     n.genes.sample = 1000, 
                                     n.iter = 1,
                                     param.list = list(A = "A", B = "B"), 
                                     n.cores = 4, 
                                     lamian.iter = 500)
View(lamian_res$Lamian_results_tidy[[1]])
# AUC-ROC
yardstick::roc_auc(lamian_res$Lamian_results_tidy[[1]], 
                   truth = factor(geneStatus_overall, levels = c("Dynamic", "NotDynamic")), 
                   .estimate = pvalue, 
                   event_level = "second")
roc_curve <- yardstick::roc_curve(lamian_res$Lamian_results_tidy[[1]], 
                                  truth = factor(geneStatus_overall, levels = c("Dynamic", "NotDynamic")), 
                                  .estimate = pvalue, 
                                  event_level = "second")
ggplot(roc_curve, aes(x = 1 - specificity, y = sensitivity)) + 
  geom_line(size = 1, alpha = 0.8) + 
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed", size = 1) + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  labs(x = "1 - Specificity", 
       y = "Sensitivity", 
       title = "ROC-AUC") + 
  theme_classic(base_size = 15)
# accuracy
yardstick::accuracy(lamian_res$Lamian_results_tidy[[1]], 
                    truth = factor(geneStatus_overall, levels = c("Dynamic", "NotDynamic")), 
                    estimate = factor(ifelse(gene_dynamic_overall == 1, "Dynamic", "NotDynamic"), 
                                      levels = c("Dynamic", "NotDynamic")))
# NIR 
mean(lamian_res$Lamian_results_tidy[[1]]$geneStatus_overall == "NotDynamic")


cell_anno <- colData(lung_sim_DEG_20_CELLS_100_balanced) %>% 
             as.data.frame() %>% 
             mutate(Cell = rownames(.), .before = 1) %>% 
             select(Cell, Sample = subject) %>% 
             mutate(across(everything(), as.character)) %>% 
             magrittr::set_rownames(NULL)
cell_pt <- as.integer(rank(lung_sim_DEG_20_CELLS_100_balanced$cell_time_normed))
names(cell_pt) <- cell_anno$Cell
samp_design <- data.frame(intercept = rep(1, length(unique(cell_anno$Sample)))) %>%
               magrittr::set_rownames(unique(cell_anno$Sample)) %>%
               as.matrix()
cell_expr <- as.matrix(logcounts(lung_sim_DEG_20_CELLS_100_balanced))[sample(1:nrow(lung_sim_DEG_20_CELLS_100_balanced), 10), ]
lamian_res <- lamian_test(expr = cell_expr,
                          cellanno = cell_anno,
                          pseudotime = cell_pt,
                          design = samp_design,
                          test.type = "time", 
                          test.method = "permutation", 
                          permuiter = 20, 
                          ncores = 2)
View(lamian_res$statistics)

global_test_results <- lamian_res$statistics %>% 
                       arrange(pval.overall) %>% 
                       mutate(pvalue_adj = stats::p.adjust(pval.overall, method = "bonferroni")) %>% 
                       dplyr::rename(pvalue = pval.overall, 
                                     Zstat = z.overall) %>% 
                       dplyr::mutate(gene = rownames(.), 
                                     gene_dynamic_overall = dplyr::case_when(pvalue_adj < 0.01 ~ 1, TRUE ~ 0)) %>% 
                       dplyr::select(gene, 
                                     Zstat, 
                                     pvalue,
                                     pvalue_adj, 
                                     gene_dynamic_overall) %>% 
                       dplyr::inner_join((rowData(lung_sim_DEG_20_CELLS_100_balanced) %>%
                                          as.data.frame() %>%
                                          dplyr::select(geneStatus_overall) %>% 
                                          dplyr::mutate(gene = rownames(.))), 
                                         by = "gene")
View(global_test_results)

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

lamian.res <- lamian_res

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
  populationFit <- t(get_pop_fit_Lamian(lamian.res = lamian_res) %>% purrr::reduce(cbind))
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

resp_preds <- compute_raw_preds_Lamian(lamian.preds = lamian_preds, orig.sce = lung_sim_DEG_20_CELLS_100_balanced)

# testing w/ raw & normalized counts 
true_vals_norm <- cell_expr["GABRE", , drop = FALSE]
true_vals_raw <- counts(lung_sim_DEG_20_CELLS_100_balanced)["GABRE", , drop = FALSE]
size_factors <- sizeFactors(lung_sim_DEG_20_CELLS_100_balanced)

est_vals_norm <- log2(t(t(true_vals_raw) / size_factors) + 1)
cor(est_vals_norm[1, ], true_vals_norm[1, ])

est_vals_raw <- t(t(2^(est_vals_norm) - 1) * size_factors)
cor(est_vals_raw[1, ], true_vals_raw[1, ])

# testing on predicted values 
lamian_preds <- get_preds_Lamian(lamian.res = lamian_res)
rownames(lamian_preds)

true_vals_raw <- counts(lung_sim_DEG_20_CELLS_100_balanced)["FUT11", , drop = FALSE]
true_vals_norm <- cell_expr["FUT11", , drop = FALSE]
pred_vals <- lamian_preds["FUT11", , drop = FALSE]
cor(pred_vals[1, ], true_vals_raw[1, ])
plot(pred_vals[1, ], true_vals_raw[1, ])

inv_true_vals <- log2(t(t(true_vals_raw) / size_factors) + 1)
cor(inv_true_vals[1, ], true_vals_norm[1, ])
cor(inv_true_vals[1, ], true_vals_raw[1, ])
plot(inv_true_vals[1, ], true_vals_raw[1, ])
cor(inv_true_vals[1, ], pred_vals[1, ])

inv_pred_vals <- t(t(2^(pred_vals) - 1) * size_factors)
cor(inv_pred_vals[1, ], true_vals_raw[1, ])
plot(inv_pred_vals[1, ], true_vals_raw[1, ])


yardstick::rmse_vec(truth = counts(lung_sim_DEG_20_CELLS_100_balanced)["GABRE", ], estimate = resp_preds["GABRE", ])
