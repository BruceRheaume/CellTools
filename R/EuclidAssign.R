#' Determine Euclidean Distance
#'
#' This function computes distances in multidimensional space and assigns variables based nearest neighbors
#'
#' @param ref Reference UMAP
#' @param query Query UMAP
#' @param label The reference label to transfer to query cells
#' @return Returns a data frame containing the distance to the nearest reference point for each query point and the label of that reference cell.
#' @examples
#'# Create random test reference
#'test.ref <- data.frame(UMAP_1 = c(runif(5, min= -5, max = -3), runif(5, min= 5, max = 8), runif(5, min= 3, max = 6)),
#'                       UMAP_2 = c(runif(5, min = -5, max = -3), runif(5, min = 5, max = 6), runif(5, min= -1, max = 3))) 
#'# Create random labels for reference
#'Col <- as.factor(c(rep(1, 5), rep(2, 5), rep(3, 5))) 
#'
#'# Create random test query
#'test.query <- data.frame(UMAP_1 = runif(5, min= -10, max = 10), UMAP_2 = runif(5, min= -10, max = 10)) 
#'
#'# Plot cells with labels to view layout
#'pal <- rainbow(3)
#'plot(test.query, pch = 16, col = "black")
#'points(test.ref, pch = 16, col = pal[Col])
#'
#'# Run EuclidAssign
#'EuclidAssign(query = test.query, ref = test.ref, label = Col)
#'
#' @references  
#' Rheaume, B. A., & Trakhtenberg, E. F. (2021). Self-learning optimization algorithm for integration 
#' of scRNA-seq datasets improves the identification of resilient and susceptible retinal ganglion cell types. 
#' \href{https://www.biorxiv.org/content/10.1101/2021.10.15.464552v2.abstract}{bioRxiv}.
#' @export
EuclidAssign <- function(query, ref, label){
  dist <- function(xy1,xy2 = ref,lab = label){
    a <- as.numeric(apply(xy2[,c(1,2)], 1, function(x){
      return(sqrt(((abs(x[1]-as.numeric(xy1)[1]))^2) + ((abs(x[2]-as.numeric(xy1)[2]))^2)))
    }))
    return(data.frame(Ref_cell = rownames(xy2)[which(a == min(a))],
                      Label = label[which(a == min(a))],
                      Dist = min(a)))
  }
  dat <- do.call(rbind,apply(query, 1, dist))
  lines <- rbind(ref[as.numeric(dat$Ref_cell),], query)
  lines$group <- rep(seq(1:nrow(query)), 2)
  print(ggplot2::ggplot(test.ref) + ggplot2::geom_point(aes(x = UMAP_1, y = UMAP_2, col = label)) + ggplot2::theme_classic() +
    ggplot2::geom_point(data = test.query, aes(x = UMAP_1, y = UMAP_2)) +
    ggplot2::geom_line(data = lines, aes(x = UMAP_1, y = UMAP_2, group = group, col = rep(dat$Label, 2)), linetype = "dashed"))
  return(dat)
}




