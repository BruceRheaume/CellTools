#' Map query Seurat object to reference Seurat object
#'
#' This function is a wrapper for Seurat's TransferData, IntegrateEmbeddings, and ProjectUMAP functions, meant to be run after
#' OptiK, as part of the Denoise/OptiK/MapTo pipeline.
#'
#' @param ref Reference Seurat object 
#' @param query Wuery Seurat object
#' @param lab The reference label to treat as clusters. Default is "seurat_clusters"
#' @param k The k-weight to use when weighting anchors. Default is opti
#' @param get.anchors A logical that when true will first run Seurat's FindTransferAnchors. Default is TRUE.
#' @param dims Number of PCs to use when get.anchors = T. Default is 30
#' @return Returns a Seurat object with a new slot where the k-weight is stored: misc$CellTools$opti_k
#' @references  
#' Rheaume, B. A., & Trakhtenberg, E. F. (2022). Self-learning algorithm for denoising 
#' and advancing the integration of scRNA-seq datasets improves the identification of resilient 
#' and susceptible retinal ganglion cell types 
#' \href{https://www.biorxiv.org/content/10.1101/2021.10.15.464552v2.abstract}{bioRxiv}.
#' @export 
MapTo <- function(ref, query, lab = "seurat_clusters", k = "opti", get.anchors = T, dims = 20, anchors = NA, neighbors = 30, nfeatures = 2000){  
  
  if(k == "opti"){
    k <- ref@misc$CellTools$opti_k
  } else {k <- k}
  
  if(get.anchors == T){
    features <- Seurat::SelectIntegrationFeatures(object.list = list(ref,query), nfeatures = nfeatures)
    anchors <- Seurat::FindTransferAnchors(
      reference = ref,
      query = query,
      features = features,
      normalization.method = "LogNormalize",
      reference.reduction = "pca",
      reduction = "pcaproject",
      dims = 1:dims
  )
  } else {anchors <- anchors}
  dat <- Seurat::TransferData(
    anchorset = anchors,
    reference = ref,
    query = query,
    refdata = lab,
    k.weight = k,
    weight.reduction = "pcaproject",
    verbose = F
  )
  dat <- Seurat::IntegrateEmbeddings(
    anchorset = anchors,
    reference = ref,
    query = dat, 
    new.reduction.name = "refpca",
    reductions = "pcaproject",
    k.weight = k,
    verbose = F
  )
  dat <- Seurat::ProjectUMAP(
    query = dat, 
    query.reduction = "refpca", 
    reference = ref, 
    reference.reduction = "pca", 
    reduction.model = "umap",
    #k.weight = k,
    #n.trees = n,
    verbose = F,
    k.param = neighbors
  )
  return(dat)
}
