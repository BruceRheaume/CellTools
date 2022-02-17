#' Optimization of threshold for capturing bad cells
#'
#' This function takes query and reference Seurat objects as input, identifies the optimal k-weight-parameter, 
#' and outputs the query object with updated metadata and miscellaneous slots, based on arguments; see below.
#'
#' @param ref Reference Seurat object
#' @param query Query Seurat object
#' @param lab The reference @meta.data slot to transfer to query object
#' @param range The range of k-weight parameters to loop through; default: 10-100
#' @param ref.dims The number of PCs to use when creating the reference UMAP and discovering anchors.
#' @param query.dims The number of PCs to use when creating the query UMAP 
#' @return A Seurat object containing the transferred data labels stored in "predicted.ids," the 
#' confidence scores saved in predicted.id.score, and the reference umap of the query object.
#' @references  
#' Rheaume, B. A., & Trakhtenberg, E. F. (2021). Self-learning optimization algorithm for integration 
#' of scRNA-seq datasets improves the identification of resilient and susceptible retinal ganglion cell types. 
#' \href{https://www.biorxiv.org/content/10.1101/2021.10.15.464552v2.abstract}{bioRxiv}.
#' @export
RefOptiThresh <- function(ref, query, lab, range = c(10,100), ref.dims = 30, query.dims = 30){
  opti_k <- ref@misc$CellTools$opti_k
  
  # Split again, but perhaps store the split values in RefOptiK?
  ref$split <- 1
  set.seed(1984)
  ref$split[sample(1:ncol(ref), ncol(query), replace = F)] <- 2
  ref.list <- Seurat::SplitObject(ref, split.by = "split")
  
  predicted.bad.df <<- data.frame(score = ref.list[[2]]$predicted.id.score, correct = rep("correct", ncol(ref.list[[2]])))
  predicted.bad.df$correct[which(ref.list[[2]]$predicted.id != ref.list[[2]]$type)] <- "incorrect"
  predicted.bad.df$correct <- as.factor(predicted.bad.df$correct)
  
  
  
  return(ref)
}