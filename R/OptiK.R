#' Find the optimal k-weight parameter for KNN 
#'
#' This function iterates through a select range of k-weight parameters and stores the optimal parameter in the Seurat object
#'
#' @param data Reference Seurat object 
#' @param lab The reference label to treat as clusters. Default is "seurat_clusters"
#' @param range The range of k-weight values to loop through. Default is 10-100.
#' @param dims The number of PCs to consider when discovering anchors
#' @param perc The percent of cells to subsample from the refeerence for the test query. Default is 20%.
#' @param num If perc = F, OptiK will look for a specific number of cells specified here.
#' @param seed Set the seed for reproducible results. Default = 1984.
#' @return Returns a Seurat object with a new slot where the k-weight is stored: misc$CellTools$opti_k
#' @references  
#' Rheaume, B. A., & Trakhtenberg, E. F. (2022). Self-learning algorithm for denoising 
#' and advancing the integration of scRNA-seq datasets improves the identification of resilient 
#' and susceptible retinal ganglion cell types 
#' \href{https://www.biorxiv.org/content/10.1101/2021.10.15.464552v2.abstract}{bioRxiv}.
#' @export 
OptiK <- function(data, lab = "seurat_clusters", range = c(10,100), dims = 30, perc = .20, num = NA, seed = 1984, n = 30){
  data@misc$CellTools <- list() 
  
  if(perc == F){
    num <- num
  } else {num <- round(ncol(data)*perc)}
  
  data$split <- 1
  set.seed(seed)
  data$split[sample(1:ncol(data), num, replace = F)] <- 2
  data.list <- Seurat::SplitObject(data, split.by = "split")
  
  # Get UMAP model
  data.list[[1]] <- Seurat::RunPCA(data.list[[1]], npcs = dims)
  data.list[[1]] <- Seurat::RunUMAP(data.list[[1]], dims = 1:dims, return.model = T)
  
  features <- Seurat::SelectIntegrationFeatures(object.list = data.list, nfeatures = 2000)
  
  anchors <- Seurat::FindTransferAnchors(
    reference = data.list[[1]],
    query = data.list[[2]],
    features = features,
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    reduction = "pcaproject",
    dims = 1:dims
  )
  optimize_k <- data.frame(k = 0, correct = 0)
  pb <- txtProgressBar(min = range[1],     
                       max = range[2], 
                       style = 3,    
                       width = 100,   
                       char = "=")
  for(k in range[1]:range[2]){ 
    data.list[[2]] <- MapTo(ref = data.list[[1]], query = data.list[[2]], lab = data.list[[1]][[lab]][,1], k = k, get.anchors = F, dims = NA, anchors = anchors, nfeatures = 2000, neighbors = n)
    optimize_new <- data.frame(k = k, correct = length(which(data.list[[2]]$predicted.id == data.list[[2]][[lab]][,1]))/ncol(data.list[[2]]))
    optimize_k <- rbind(optimize_new, optimize_k)
    setTxtProgressBar(pb, k)
  }
  close(pb)
  optimize_k <- optimize_k[-nrow(optimize_k),]
  opti_k <- min(optimize_k$k[which(optimize_k$correct == max(optimize_k$correct))]) #18
  #optimized_k <- ggplot2::ggplot(optimize_k) + ggplot2::geom_line(ggplot2::aes(x = k, y = correct)) + ggplot2::geom_vline(xintercept = opti_k, col = "red")
  #print(optimized_k)
  
  #data.list[[2]] <- MapTo(ref = data.list[[1]], query = data.list[[2]], label = data.list[[1]][[lab]][,1], k = opti_k)
  
  # Save the data in Seurat object slot
  # data$mapto.predicted <- 0
  # data$mapto.predicted[which(data$split == 2)] <- data.list[[2]]$predicted.id
  # data$mapto.scores <- 0
  # data$mapto.scores[which(data$split == 2)] <- data.list[[2]]$predicted.id.score
  
  data@misc$CellTools$opti_k <- opti_k
  data@misc$CellTools$optimize_k <- optimize_k
  #data@misc$CellTools$optimized_k <- optimized_k
  
  #data@misc$CellTools$ref.umap <- data.list[[1]]@reductions$umap@cell.embeddings
  #data@misc$CellTools$query.umap <- data.list[[2]]@reductions$ref.umap@cell.embeddings
  
  return(data)
}
