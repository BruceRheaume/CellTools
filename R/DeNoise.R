#' Denoise a reference Seurat Object
#'
#' This function denoises, or removes sources of noise, which can affect downstream intergration analyses
#'
#' @param ref Reference Seurat object with pre-calculated UMAP
#' @param lab The reference label to treat as clusters. Default is "seurat_clusters"
#' @param coef The coefficient to multiply by the IQR. Default is 1.5
#' @return Returns a Seurat object with a new slot called "denoise"
#' @references  
#' Rheaume, B. A., & Trakhtenberg, E. F. (2022). Self-learning algorithm for denoising 
#' and advancing the integration of scRNA-seq datasets improves the identification of resilient 
#' and susceptible retinal ganglion cell types 
#' \href{https://www.biorxiv.org/content/10.1101/2021.10.15.464552v2.abstract}{bioRxiv}.
#' @export 
DeNoise <- function(ref, lab, coef = 1.5){
  ref.embed <- as.data.frame(ref@reductions$umap@cell.embeddings)
  ref.embed$cluster <- ref[[lab]][,1]
  
  cls.x.l <- lapply(unique(ref.embed$cluster), function(x){
    cls.x <- ref.embed$UMAP_1[which(ref.embed$cluster == x)]
    return(data.frame(cluster = rep(x, length(boxplot.stats(cls.x, coef = coef)$out)),
                      value = boxplot.stats(cls.x, coef = coef)$out))
  })
  cls.x <- do.call(rbind, cls.x.l)
  
  cois <- mapply(function(cluster, value){
    return(rownames(ref.embed)[which(ref.embed$UMAP_1 == value & ref.embed$cluster == cluster)])
  },cluster = cls.x$cluster, value = cls.x$value)
  
  ref.embed2 <- ref.embed[-which(rownames(ref.embed) %in% cois),]
  
  cls.y.l <- lapply(unique(ref.embed2$cluster), function(x){
    cls.y <- ref.embed2$UMAP_2[which(ref.embed2$cluster == x)]
    return(data.frame(cluster = rep(x, length(boxplot.stats(cls.y, coef = coef)$out)),
                      value = boxplot.stats(cls.y, coef = coef)$out))
  })
  cls.y <- do.call(rbind, cls.y.l)
  
  cois <- mapply(function(cluster, value){
    return(rownames(ref.embed2)[which(ref.embed2$UMAP_2 == value & ref.embed2$cluster == cluster)])
  },cluster = cls.y$cluster, value = cls.y$value)
  
  ref.embed3 <- ref.embed2[-which(rownames(ref.embed2) %in% cois),]
  
  print(ggplot2::ggplot(ref.embed3) + ggplot2::geom_point(ggplot2::aes(x = UMAP_1, y = UMAP_2, col = cluster), size = 0.5) + 
          ggplot2::theme_classic())
  
  ref$denoise <- "noise"
  ref$denoise[which(colnames(ref) %in% rownames(ref.embed3))] <- "signal"
  ref$denoise <- as.factor(ref$denoise)
  return(ref)
}




