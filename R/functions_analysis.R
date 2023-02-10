# downstream analysis functions

create_metric_table <- function(metric.list = NULL) {
  # check inputs 
  if (is.null(metric.list)) { stop("you must supply a list of model metrics to create_metric_table().") }
  n_cols <- purrr::map_dbl(metric.list, ncol)
  if (length(unique(n_cols)) > 1) { stop("metric tables must all have the same dimensions.") }
  
  # collate tables
  master_metric_table <- purrr::reduce(metric.list, rbind) %>% 
                         dplyr::mutate(MASTER_DATE_GENERATED = Sys.Date())
  return(master_metric_table)
}

save_figures <- function(fig.object = NULL, 
                         fig.name = NULL, 
                         fig.dir = NULL, 
                         fig.dim = c(8, 6), 
                         fig.device = "png", 
                         fig.units = "in", 
                         fig.dpi = "retina") {
  # check inputs 
  if (is.null(fig.object) || is.null(fig.name) || is.null(fig.dir)) { stop("fig.object, fig.name, & fig.dir must be non-NULL") }
  fig.dir.full <- paste0("/blue/rbacher/j.leary/repos/scLANE-Sims/Figures/", fig.dir, "/")
  if (!dir.exists(fig.dir.full)) {
    dir.create(fid.dir.full)
  }
  fig.name.full <- paste0(fig.name, ".", fig.device)
  
  # save ggplot2 figure
  ggplot2::ggsave(filename = fig.name.full,
                  plot = fig.object, 
                  device = fig.device, 
                  path = fig.dir.full, 
                  width = fig.dim[1],
                  height = fig.dim[2], 
                  units = fig.units, 
                  dpi = fig.dpi)
}
