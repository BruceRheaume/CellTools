'%!in%' <- function(x,y)!('%in%'(x,y))

#' Optimized Tracing  of Query to Reference Data
#'
#' Stupid Problem Solver
#'
#' @param ref Reference Seurat object
#' @param query Query Seurat object
#' @param lab The reference @meta.data slot to transfer to query object
#' @param method Whether or not to do the reintegration. If method = "one-step," only the k-parameter optimization and standard mapping will be performed (MapTo). If method = "rescue," the confidence-score optimization and reintegration will also be performed; this step greatly increases the accuracy of predicted cell IDs (MapToRescue); default is "rescue."
#' @param range The range of k-weight parameters to loop through; default: 10-100
#' @param sample If a number is given, the reference and query data will be downsampled to the sample size given by this parameter. This will vastly increase the speed of the algorithm, but it will decrease the accuracy. Default is FALSE.
#' @param dims The number of PCs to use when creating the reference UMAP and when discovering anchors.
#' @return A Seurat object containing the transferred data labels stored in "predicted.ids," the 
#' confidence scores saved in predicted.id.score, and the reference umap of the query object.
#' @references  
#' Rheaume, B. A., & Trakhtenberg, E. F. (2021). Self-learning optimization algorithm for integration 
#' of scRNA-seq datasets improves the identification of resilient and susceptible retinal ganglion cell types. 
#' \href{https://www.biorxiv.org/content/10.1101/2021.10.15.464552v2.abstract}{bioRxiv}.
#' @export
map_to <- function(ref, query, lab, method = "Rescue", range = c(10,100), sample = F, dims = 30){
  map.to <- function(ref = ref, query = query, label = lab, k = 30){
    n = 200  
    dat <- Seurat::TransferData(
      anchorset = anchors,
      reference = ref,
      query = query,
      refdata = label,
      k.weight = k,
      weight.reduction = "pcaproject"
    )
    dat <- Seurat::IntegrateEmbeddings(
      anchorset = anchors,
      reference = ref,
      query = dat, 
      new.reduction.name = "ref.pca",
      reductions = "pcaproject",
      k.weight = k
    )
    dat <- Seurat::ProjectUMAP(
      query = dat, 
      query.reduction = "ref.pca", 
      reference = ref, 
      reference.reduction = "pca", 
      reduction.model = "umap",
      k.weight = k,
      n.trees = n
    )
    return(dat)
  }
  
  if(sample != F){
    set.seed(1984)
    sampled.cells <- sample(x = colnames(query), size = round((sample/ncol(ref))*ncol(query)), replace = F) # query should be smaller than  reference
    query <- subset(x = query, cells = sampled.cells)
    
    set.seed(1984)
    sampled.cells <- sample(x = colnames(ref), size = sample, replace = F)
    ref <- subset(x = ref, cells = sampled.cells)
  } 
  
  # First split the reference to test query and test reference
  ref$split <- 1
  set.seed(1984)
  ref$split[sample(1:ncol(ref), ncol(query), replace = F)] <- 2
  ref.list <- Seurat::SplitObject(ref, split.by = "split")
  features <- Seurat::SelectIntegrationFeatures(object.list = ref.list, nfeatures = 2000)
  ref.list <- lapply(ref.list, Seurat::RunPCA)
  ref.list[[1]] <- Seurat::RunUMAP(ref.list[[1]], dims = 1:dims, return.model = T) 
  return(ref.list)
}
