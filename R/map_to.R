#' Map Query Data to Reference UMAP
#'
#' This function takes Seurat query and reference objects and identifies the optimal k-weight-parameter 
#' and outputs the query object with updated metadata; see below.
#'
#' @param ref Reference Seurat object
#' @param query Query Seurat object
#' @param method The method to use for k-weight optimization. Options include: "id.score" (default), "sub-sample"
#' @param range The range of k-weight parameters to loop through; default: 5-10
#' @return A Seurat object containing the transferred data labels stored in "predicted.ids," the 
#' confidence scores saved in predicted.id.score, and the reference umap of the query object.
#' @references  
#' Rheaume, B. A., & Trakhtenberg, E. F. (2021). Self-learning optimization algorithm for integration 
#' of scRNA-seq datasets improves the identification of resilient and susceptible retinal ganglion cell types. 
#' \href{https://www.biorxiv.org/content/10.1101/2021.10.15.464552v2.abstract}{bioRxiv}.
#' @export
map_to <- function(ref, query, method = "id.score", range = c(5,10)){
  in.dt <- data.table::fread(infile, header = TRUE)
  in.dt <- in.dt[!duplicated(in.dt[, 1]), ]
  in.mat <- as.matrix(in.dt[, -1, with = FALSE])
  rownames(in.mat) <- unlist(in.dt[, 1, with = FALSE])
  in.mat
}