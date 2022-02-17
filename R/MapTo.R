#' Optimization of K-Nearest Neighbors for Tracing of Query to Reference Data
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
MapTo <- function(ref, query, lab, range = c(10,100), ref.dims = 30, query.dims = 30){
  query@misc$CellTools <- list() 
  # MapTo:
  map.to <- function(ref = ref, query = query, label = lab, k = 30){
    n = 200  
    dat <- Seurat::TransferData(
      anchorset = anchors,
      reference = ref,
      query = query,
      refdata = label,
      k.weight = k,
      weight.reduction = "pcaproject",
      verbose = F
    )
    dat <- Seurat::IntegrateEmbeddings(
      anchorset = anchors,
      reference = ref,
      query = dat, 
      new.reduction.name = "refpca_",
      reductions = "pcaproject",
      k.weight = k,
      verbose = F
    )
    dat <- Seurat::ProjectUMAP(
      query = dat, 
      query.reduction = "refpca_", 
      reference = ref, 
      reference.reduction = "pca", 
      reduction.model = "umap",
      #k.weight = k,
      n.trees = n,
      verbose = F
    )
    return(dat)
  }
  
  # 1st step: split the reference
  ref$split <- 1
  set.seed(1984)
  ref$split[sample(1:ncol(ref), ncol(query), replace = F)] <- 2
  ref.list <- Seurat::SplitObject(ref, split.by = "split")
  
  # Get UMAP model
  #ref.list[[2]] <- Seurat::RunPCA(ref.list[[2]], npcs = ref.dims)
  features <- Seurat::SelectIntegrationFeatures(object.list = ref.list, nfeatures = 2000)
  
  anchors <- Seurat::FindTransferAnchors(
    reference = ref.list[[1]],
    query = ref.list[[2]],
    features = features,
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    reduction = "pcaproject",
    dims = 1:ref.dims
  )
  optimize_k <- data.frame(k = 0, correct = 0)
  pb <- txtProgressBar(min = range[1],      # Minimum value of the progress bar
                       max = range[2], # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")
  for(k in range[1]:range[2]){ 
    ref.list[[2]] <- map.to(ref = ref.list[[1]], query = ref.list[[2]], label = ref.list[[1]][[lab]][,1], k = k)
    optimize_new <- data.frame(k = k, correct = length(which(ref.list[[2]]$predicted.id == ref.list[[2]][[lab]][,1]))/ncol(ref.list[[2]]))
    optimize_k <- rbind(optimize_new, optimize_k)
    setTxtProgressBar(pb, k)
  }
  close(pb)
  optimize_k <- optimize_k[-nrow(optimize_k),]
  opti_k <<- min(optimize_k$k[which(optimize_k$correct == max(optimize_k$correct))]) #18
  optimized_k <- ggplot2::ggplot(optimize_k) + ggplot2::geom_line(ggplot2::aes(x = k, y = correct)) + ggplot2::geom_vline(xintercept = opti_k, col = "red")
  print(optimized_k)
  
  # Save the data in Seurat object slot
  query@misc$CellTools$opti_k <- opti_k
  query@misc$CellTools$optimize_k <- optimize_k
  query@misc$CellTools$optimized_k <- optimized_k
  
  ref.list[[2]] <- map.to(ref = ref.list[[1]], query = ref.list[[2]], label = ref.list[[1]][[lab]][,1], k = opti_k)
  
  return(query)
}