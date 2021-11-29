'%!in%' <- function(x,y)!('%in%'(x,y))

#' Map Query Data to Reference UMAP
#'
#' This function takes Seurat query and reference objects and identifies the optimal k-weight-parameter 
#' and outputs the query object with updated metadata; see below.
#'
#' @param ref Reference Seurat object
#' @param query Query Seurat object
#' @param lab The reference meta.data slot to transfer to query object
#' @param method The method to use for k-weight optimization. Options include: "id.score" (default) or "ref-sample"
#' @param range The range of k-weight parameters to loop through; default: 5-100
#' @param sample If method = "id.score," the reference and query data will be downsampled to the sample size given by this parameter. If method = "ref-sample," the sample number will be the size of the "query" object split from the reference.
#' @param dims Number of PCs to use when creating the reference UMAP and when discovering anchors.
#' @return A Seurat object containing the transferred data labels stored in "predicted.ids," the 
#' confidence scores saved in predicted.id.score, and the reference umap of the query object.
#' @references  
#' Rheaume, B. A., & Trakhtenberg, E. F. (2021). Self-learning optimization algorithm for integration 
#' of scRNA-seq datasets improves the identification of resilient and susceptible retinal ganglion cell types. 
#' \href{https://www.biorxiv.org/content/10.1101/2021.10.15.464552v2.abstract}{bioRxiv}.
#' @export
map_to <- function(ref, query, lab, method = "id.score", range = c(5,100), sample = F, dims = 50){
  if(method == "id.score"){
    if(sample != F){
      set.seed(1984)
      sampled.cells <- sample(x = colnames(ref), size = sample, replace = F)
      ref2 <- subset(x = ref, cells = sampled.cells)
      
      set.seed(1984)
      sampled.cells <- sample(x = colnames(query), size = sample, replace = F)
      query2 <- subset(x = query, cells = sampled.cells)
      
      #lab2 <- lab[which(names(lab) %in% colnames(ref2))]
    } else {ref2 <- ref; query2 <- query; } #lab2 <- lab
    set.seed(1984)
    ref2 <- Seurat::RunPCA(ref2)
    ref2 <- Seurat::RunUMAP(ref2, reduction = "pca", dims = 1:dims, return.model = T)
    
    features <- Seurat::SelectIntegrationFeatures(object.list = list(ref2, query2), nfeatures = 2000)
    
    anchors <- Seurat::FindTransferAnchors(
      reference = ref2,
      query = query2,
      features = features,
      normalization.method = "LogNormalize",
      reference.reduction = "pca",
      reduction = "rpca",
      dims = 1:dims
    )
    iterate.list <- list()
    for(k in range[1]:range[2]){
      query2 <- Seurat::TransferData(
        anchorset = anchors,
        reference = ref2,
        query = query2,
        refdata = ref2[[lab]][,1],
        k.weight = k,
        weight.reduction = "pcaproject"
      )
      query2 <- Seurat::IntegrateEmbeddings(
        anchorset = anchors,
        reference = ref2,
        query = query2, 
        new.reduction.name = "ref.pca",
        reductions = "pcaproject",
        k.weight = k
      )
      query2 <- Seurat::ProjectUMAP(
        query = query2, 
        query.reduction = "ref.pca", 
        reference = ref2, 
        reference.reduction = "pca", 
        reduction.model = "umap",
        k.weight = k
      )
      iterate.list[[paste(k)]] <- data.frame(Cluster = query2$predicted.id, Score = query2$predicted.id.score)
    }
    
    score.list <- lapply(iterate.list, function(x){mean(x$Score)})
    names(score.list) <- names(iterate.list)
    score.list <- unlist(score.list)
    print(paste("Optimal k = ", names(score.list)[which(score.list == max(score.list))],
                "; Mean score = ", max(score.list), sep = ""))
    print("Running map_to on query data...")
    
    k <- as.numeric(names(score.list)[which(score.list == max(score.list))])
    
    score.df <- as.data.frame(score.list)
    score.df$k <- factor(rownames(score.df), levels = rownames(score.df))
    query@misc <- score.df
    ggplot2::ggplot(score.df, aes(x = k, y = score.list, group = 1)) + geom_line()
    
  } else if(method == "ref-sample"){
    if(sample %!in% c(50:10000)){ stop("Must enter sample size between 50 and 10,000 when method = \'ref-sample\'")
    } else {
      ref$split <- 1
      set.seed(1984)
      ref$split[sample(1:ncol(ref), sample, replace = F)] <- 2
      
      ref.list <- Seurat::SplitObject(ref, split.by = "split")
      features <- Seurat::SelectIntegrationFeatures(object.list = ref.list, nfeatures = 2000)
      ref.list <- lapply(ref.list, Seurat::RunPCA)
      ref.list[[1]] <- Seurat::RunUMAP(ref.list[[1]], dims = 1:50, return.model = T) 
      #lab2 <- lab[which(names(lab) %in% colnames(ref.list[[1]]))]
      
      anchors <- Seurat::FindTransferAnchors(
        reference = ref.list[[1]],
        query = ref.list[[2]],
        features = features,
        normalization.method = "LogNormalize",
        reference.reduction = "pca",
        reduction = "rpca",
        dims = 1:dims
      )
      scores.list <- data.frame(k = 0, correct = 0, score = 0)
      for(k in range[1]:range[2]){
        ref.list[[2]] <- Seurat::TransferData(
          anchorset = anchors,
          reference = ref.list[[1]],
          query = ref.list[[2]],
          refdata = ref.list[[1]][[lab]][,1],
          k.weight = k,
          weight.reduction = "pcaproject"
        )
        ref.list[[2]] <- Seurat::IntegrateEmbeddings(
          anchorset = anchors,
          reference = ref.list[[1]],
          query = ref.list[[2]], 
          new.reduction.name = "ref.pca",
          reductions = "pcaproject",
          k.weight = k
        )
        ref.list[[2]] <- Seurat::ProjectUMAP(
          query = ref.list[[2]], 
          query.reduction = "ref.pca", 
          reference = ref.list[[1]], 
          reference.reduction = "pca", 
          reduction.model = "umap",
          k.weight = k
        )
        scores <- c(k, length(which(ref.list[[2]]$predicted.id == ref.list[[2]][[lab]][,1]))/length(ref.list[[2]][[lab]][,1]), mean(ref.list[[2]]$predicted.id.score))
        scores.list <- rbind(scores.list, scores)
      }
    }
    print(paste("Optimal k = ", scores.list$k[which(scores.list$correct == max(scores.list$correct))],
                "; Mean score = ", scores.list$score[which(scores.list$correct == max(scores.list$correct))], sep = ""))
    print("Running map_to on query data...")
    
    k <- scores.list$k[which(scores.list$correct == max(scores.list$correct))]

    query@misc <- scores.list
    ggplot2::ggplot(scores.list, aes(x = k, y = correct, group = 1)) + geom_line()
  }
  
  features <- Seurat::SelectIntegrationFeatures(object.list = list(ref, query), nfeatures = 2000)
  anchors <- Seurat::FindTransferAnchors(
    reference = ref,
    query = query,
    features = features,
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    reduction = "rpca",
    dims = 1:dims
  )
  query <- Seurat::TransferData(
    anchorset = anchors,
    reference = ref,
    query = query,
    refdata = ref[[lab]][,1],
    k.weight = k,
    weight.reduction = "pcaproject"
  )
  query <- Seurat::IntegrateEmbeddings(
    anchorset = anchors,
    reference = ref,
    query = query, 
    new.reduction.name = "ref.pca",
    reductions = "pcaproject",
    k.weight = k
  )
  query <- Seurat::ProjectUMAP(
    query = query, 
    query.reduction = "ref.pca", 
    reference = ref, 
    reference.reduction = "pca", 
    reduction.model = "umap",
    k.weight = k
  )
  return(query)
}
