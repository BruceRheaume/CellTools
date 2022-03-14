#' Optimized Tracing  of Query to Reference Data
#'
#' This function examines the global transcriptomic properties of single cell data.
#'
#' @param dat A Seurat object
#' @param lab The reference @meta.data slot on which to compare global parameters
#' @return Global correlation analyses of the clusters provided in lab parameter
#' @references  
#' Rheaume, B. A., & Trakhtenberg, E. F. (2022). Self-learning algorithm for denoising 
#' and advancing the integration of scRNA-seq datasets improves the identification of resilient 
#' and susceptible retinal ganglion cell types 
#' \href{https://www.biorxiv.org/content/10.1101/2021.10.15.464552v2.abstract}{bioRxiv}.
#' @export
GlobalProbe <- function(dat, lab){
  
  sem <- function(x){
    se <- (sd(x))/(sqrt(length(x)))
    return(se)
  }
  
  means.dat <- data.frame(matrix(0,nrow(dat@assays$RNA@data), length(unique(dat[[lab]][,1]))))
  colnames(means.dat) <- unique(dat[[lab]][,1])

  for(k in colnames(means.dat)){
    means.dat[,k] <- apply(dat@assays$RNA@data[,which(dat[[lab]][,1] == k)], 1, mean)
  }
  dat.cors <- cor(means.dat)
  mean.cors <- apply(dat.cors, 2, mean)
  dat.cors <- dat.cors[,order(mean.cors, decreasing = F)]
  
  #Plot
  df <- data.frame(Cluster = colnames(dat.cors), Cor = apply(dat.cors,2,mean), sd = apply(dat.cors, 2, sem))
  df$Cluster <- sub("_.*","",df$Cluster)
  df <- df[order(df$Cor),]
  df$Cluster <- factor(df$Cluster, levels = df$Cluster)

  p<- ggplot2::ggplot(df, aes(x=Cluster, y=Cor, position_dodge(1))) + 
    #geom_line(cex = 0.5) +
    ggplot2::geom_point()+
    ggplot2::geom_errorbar(aes(ymin=Cor-sd, ymax=Cor+sd), width = 0.5, cex = 0.5, 
                  position=position_dodge(0.001))
  p <- p + ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = .5, size = 18),
          axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 18), legend.position = "none",
          axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.ticks = element_line(), axis.title.x = element_blank()) +
    ggplot2::ylab("Mean correlation coefficient (r)") + xlab("Clusters") +
    ggplot2::coord_cartesian(expand = T) #turn off axis expansion (padding)
                    # ,ylim = c(0.89, 0.98)) #manually set limits
  print(p)
}