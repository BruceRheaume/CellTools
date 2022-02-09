'%!in%' <- function(x,y)!('%in%'(x,y))

#' Optimized Tracing of Query to Reference Data
#'
#' This function takes query and reference Seurat objects as input, identifies the optimal k-weight-parameter, 
#' and outputs the query object with updated metadata and miscellaneous slots, based on arguments; see below.
#'
#' @param ref Reference Seurat object
#' @param query Query Seurat object
#' @param lab The reference @meta.data slot to transfer to query object
#' @param method Whether or not to do the reintegration. If method = "one-step," only the k-parameter optimization and standard mapping will be performed (MapTo). If method = "rescue," the confidence-score optimization and reintegration will also be performed; this step greatly increases the accuracy of predicted cell IDs (MapToRescue); default is "rescue."
#' @param range The range of k-weight parameters to loop through; default: 10-100
#' @param sample If a number is given, the reference and query data will be downsampled to the sample size given by this parameter. This will vastly increase the speed of the algorithm, but it will decrease the accuracy. Default is FALSE.
#' @param ref.dims The number of PCs to use when creating the reference UMAP and discovering anchors.
#' @param query.dims The number of PCs to use when creating the query UMAP 
#' @return A Seurat object containing the transferred data labels stored in "predicted.ids," the 
#' confidence scores saved in predicted.id.score, and the reference umap of the query object.
#' @references  
#' Rheaume, B. A., & Trakhtenberg, E. F. (2021). Self-learning optimization algorithm for integration 
#' of scRNA-seq datasets improves the identification of resilient and susceptible retinal ganglion cell types. 
#' \href{https://www.biorxiv.org/content/10.1101/2021.10.15.464552v2.abstract}{bioRxiv}.
#' @export
MapToRescue <- function(ref, query, lab, method = "rescue", range = c(10,100), ref.dims = 30, query.dims = 30, existing.model = F){
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
  
  dist <- function(xy1,xy2=ref.umap){
    a <- as.numeric(apply(xy2, 1, function(x){
      return(sqrt(((abs(x[1]-as.numeric(xy1)[1]))^2) + ((abs(x[2]-as.numeric(xy1)[2]))^2)))
    }))
    return(data.frame(Ref_cell = rownames(ref.umap)[which(a == min(a))],
                      Cluster = ref.seu[which(a == min(a))],
                      Dist = min(a)))
  }
  
  # if(sample != F){
  #   set.seed(1984)
  #   sampled.cells <- sample(x = colnames(query), size = round((sample/ncol(ref))*ncol(query)), replace = F) # query should be smaller than  reference
  #   query <- SeuratObject::subset(x = query, cells = sampled.cells)
  # 
  #   set.seed(1984)
  #   sampled.cells <- sample(x = colnames(ref), size = sample, replace = F)
  #   ref <- subset(x = ref, cells = sampled.cells)
  # }
  
  # First split the reference to test query and test reference
  ref$split <- 1
  set.seed(1984)
  ref$split[sample(1:ncol(ref), ncol(query), replace = F)] <- 2
  ref.list <- Seurat::SplitObject(ref, split.by = "split")
  features <- Seurat::SelectIntegrationFeatures(object.list = ref.list, nfeatures = 2000)
  ref.list <- lapply(ref.list, Seurat::RunPCA)
  if(existing.model == F){
    ref.list[[1]] <- Seurat::RunUMAP(ref.list[[1]], dims = 1:ref.dims, return.model = T) 
  }
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
  optimize_k <<- optimize_k[-nrow(optimize_k),]
  opti_k <<- min(optimize_k$k[which(optimize_k$correct == max(optimize_k$correct))]) #18
  #optimized_k <<- ggplot2::ggplot(optimize_k) + ggplot2::geom_line(aes(x = k, y = correct)) + ggplot2::geom_vline(xintercept = opti_k, col = "red")
  #print(optimized_k)
  #write.csv(optimize_k, "optimize_k.csv")

#### above works ####
  if(existing.model == F){
    ref <- Seurat::RunUMAP(ref, dims = 1:ref.dims, return.model = T)
  }
  #ref.list <- lapply(list(ref,query), Seurat::RunPCA)
  features <- Seurat::SelectIntegrationFeatures(object.list = list(ref,query), nfeatures = 2000)
  anchors <- Seurat::FindTransferAnchors(
    reference = ref,
    query = query,
    features = features,
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    reduction = "pcaproject",
    dims = 1:ref.dims
  )
query <- map.to(ref = ref, query = query, label = ref[[lab]][,1], k = opti_k)
query$maptorescue <- query$predicted.id

if(ncol(query) < 500){
  return(query)
  } else{

#### Rescue ####
# set to optimal 
prop <- ncol(query)/ncol(ref)
if(prop > 0.05){
  num <- round(prop*ncol(query))
} else {num <- round(ncol(query)/4)}

Seurat::DefaultAssay(query) <- "RNA"
all.genes <- rownames(query@assays$RNA@data)
query <- Seurat::ScaleData(query, features = all.genes)
query <- Seurat::FindVariableFeatures(query)
query <- Seurat::RunPCA(query)
query <- Seurat::FindNeighbors(query, dims = 1:query.dims)
query <- Seurat::FindClusters(query, resolution = 0.5)

query$split <- 1
set.seed(1984)
query$split[sample(1:ncol(query), num, replace = F)] <- 2

query.list <- Seurat::SplitObject(query, split.by = "split")
features <- Seurat::SelectIntegrationFeatures(object.list = query.list, nfeatures = 2000)
query.list <- lapply(query.list, Seurat::RunPCA)
query.list[[1]] <- Seurat::RunUMAP(query.list[[1]], dims = 1:query.dims, return.model = T) 

anchors <- Seurat::FindTransferAnchors(
  reference = query.list[[1]],
  query = query.list[[2]],
  features = features,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  reduction = "pcaproject",
  dims = 1:query.dims
)
query.list[[2]] <- map.to(ref = query.list[[1]], query = query.list[[2]], label = query.list[[1]]$seurat_clusters, k = opti_k) 
length(which(query.list[[2]]$predicted.id != query.list[[2]]$seurat_clusters))

# Make a cutoff predicted id score for maximal incorrect capture
predicted.bad.df <<- data.frame(score = query.list[[2]]$predicted.id.score, 
                               correct = rep("correct", ncol(query.list[[2]])))
predicted.bad.df$correct[which(query.list[[2]]$predicted.id != query.list[[2]]$seurat_clusters)] <- "incorrect"
predicted.bad.df$correct <- as.factor(predicted.bad.df$correct)

#good_bad_viols <<- ggplot2::ggplot(predicted.bad.df) + ggplot2::geom_violin(aes(x = correct, y = score, fill = correct), scale = "width") +
#  ggplot2::scale_fill_manual(values = c("green", "red")) 
#print(good_bad_viols)
# make.fig(p1, "injured_correct_incorrect_predicted.ids", 5,5)


# DENSITY OVERLAP
lower.limit <- min(predicted.bad.df$score)
upper.limit <- max(predicted.bad.df$score)
correct.density <- density(subset(predicted.bad.df, correct == "correct")$score, from = lower.limit, to = upper.limit, n = 2^10)
not.correct.density <- density(subset(predicted.bad.df, correct == "incorrect")$score, from = lower.limit, to = upper.limit, n = 2^10)

density.difference <- correct.density$y - not.correct.density$y
intersection.point <<- min(correct.density$x[which(diff(density.difference > 0) != 0) + 1])

# good_bad_density <<- ggplot2::ggplot(predicted.bad.df, aes(x = score, fill = correct)) + ggplot2::geom_density(alpha = 0.7) +
#   ggplot2::scale_fill_manual(values = c("green", "red")) + theme_classic() +
#   ggplot2::geom_vline(xintercept = intersection.point, color = "red")
# print(good_bad_density)
#make.fig(p1, "correct_incorrect_densities_intersect", 5,5)

#length(which(query$predicted.id.score < intersection.point))

#### Standard integration/euclidean distance ####
bad.cells <- colnames(query)[which(query$predicted.id.score < intersection.point)]
query$split.bad <- 1
query$split.bad[which(colnames(query) %in% bad.cells)] <- 2

query.bad <- subset(query, subset = split.bad == 2)
features.bad <- Seurat::SelectIntegrationFeatures(object.list = list(ref, query.bad), nfeatures = 2000)
anchors.bad <-  Seurat::FindIntegrationAnchors(object.list = list(ref, query.bad), anchor.features = features.bad)

# this command creates an 'integrated' data assayy  
query.combined <- Seurat::IntegrateData(anchorset = anchors.bad)
query.combined <- Seurat::ScaleData(query.combined, verbose = FALSE)
query.combined <- Seurat::RunPCA(query.combined, npcs = query.dims, verbose = FALSE)
query.combined <- Seurat::RunUMAP(query.combined, reduction = "pca", dims = 1:query.dims)

new.embed <- as.data.frame(query.combined@reductions$umap@cell.embeddings)
new.embed$shape <- "reference"
new.embed[bad.cells,]$shape <- "incorrect"
new.embed$type <- query.combined[[lab]][,1]
new.embed[bad.cells,]$type <- query.bad$predicted.id

#new.embed$cluster <- query.combined$seurat_clusters
new.embed$shape <- factor(new.embed$shape, levels = c("reference", "incorrect"))
new.embed <<- new.embed

# p1 <- ggplot(new.embed) + geom_point(aes(x = UMAP_1, y = UMAP_2, col = type, shape = shape, size = shape)) +
#   scale_shape_manual(values = c(16,4)) + scale_size_manual(values = c(.5,2)) + theme_classic()
# p1
# make.fig(p1, "new_integration_bad_predict_ref_onc", 5,6)

# new_integration <<- ggplot2::ggplot(new.embed) + ggplot2::geom_point(aes(x = UMAP_1, y = UMAP_2, col = type, shape = shape, size = shape)) +
#   ggplot2::scale_shape_manual(values = c(16,4)) + scale_size_manual(values = c(.5,2)) + theme_classic()
# print(new_integration)

# Extract positions

# first break into two embeddings 
ref.umap <- query.combined@reductions$umap@cell.embeddings
ref.umap <- ref.umap[-which(rownames(ref.umap) %in% bad.cells),]
ref.seu <- query.combined[[lab]][,1][-which(colnames(query.combined) %in% bad.cells)]

bad.umap <- query.combined@reductions$umap@cell.embeddings[bad.cells,]
bad.seu <- query.combined$predicted.id[which(colnames(query.combined) %in% bad.cells)]

print(nrow(bad.umap))
#print(system.time({ apply(bad.umap[c(1:100),], 1, dist) })) # should take ~30 min


##########
new <- apply(bad.umap, 1, dist) # 30 min
new.df <- do.call(rbind, new)
new.df$Query_cell <- rownames(new.df)
new.df <- new.df[,c(1,4,2,3)]
set.seed(1984)
new.df$Rand_x <- sample(c(-1,1), nrow(new.df), replace = T)
set.seed(5)
new.df$Rand_y <- sample(c(-1,1), nrow(new.df), replace = T)

new.df$Rand_x_dist <- sapply(new.df$Dist, function(x){
  runif(1, min=0, max=x)
})

new.df$Y_comp <- mapply(function(dist, rand_x) {
  sqrt((dist^2)-(rand_x^2))
}, dist = new.df$Dist, rand_x = new.df$Rand_x_dist)

# Scale
combined.range.x <- range(query.combined@reductions$umap@cell.embeddings[,1])[2] - range(query.combined@reductions$umap@cell.embeddings[,1])[1]
ref.range.x <- range(ref@reductions$umap@cell.embeddings[,1])[2] - range(ref@reductions$umap@cell.embeddings[,1])[1]
x.scale <- ref.range.x/combined.range.x

combined.range.y <- range(query.combined@reductions$umap@cell.embeddings[,2])[2] - range(query.combined@reductions$umap@cell.embeddings[,2])[1]
ref.range.y <- range(ref@reductions$umap@cell.embeddings[,2])[2] - range(ref@reductions$umap@cell.embeddings[,2])[1]
y.scale <- ref.range.y/combined.range.y

new.df$Rand_x_dist <- new.df$Rand_x_dist*x.scale*new.df$Rand_x
new.df$Y_comp <- new.df$Y_comp*y.scale*new.df$Rand_y
new.df <<- new.df

query$predicted.id[bad.cells] <- new.df$Cluster
query$maptorescue <- quer$predicted.id

return(query)
}
}

# # Get the percentage of unassigned were "remapped"
# new.old.df2 <- new.old.df[which(a != bad.seu),]
# 
# remapped <- as.data.frame(table(new.old.df2$Type))
# remapped$orig <- sapply(remapped$Var1, function(x) length(which(query$type == x)))
# remapped$percent <- remapped$Freq/remapped$orig

