
qq <- function(...) {
  sapply(match.call()[-1], deparse)
}



plot2geneSeuratPCA = function(x,g1,g2){
  print(ggplot(data.frame(x@dr$pca@cell.embeddings), aes(x=PC1,y=PC2, color=as.numeric(x@data[g1,]), fill=as.numeric(x@data[g2,]))) + geom_point( alpha=0.7, size=2.5, shape=21, stroke=0.5) + theme(aspect.ratio = 1) + scale_colour_gradient(low="gray90",high="red") + theme_bw()  +  guides(color=guide_legend(title=g1), fill=guide_legend(title=g2)) + scale_fill_gradient(low="gray90",high="blue") + theme(aspect.ratio = 1, legend.position="bottom") + xlab("PC 1") + ylab("PC 2"))
}

plot2geneSeuratTSNE = function(x,g1,g2){
  print(ggplot(data.frame(x@dr$tsne@cell.embeddings), aes(x=tSNE_1,y=tSNE_2, color=as.numeric(x@data[g1,]), fill=as.numeric(x@data[g2,]))) +
          geom_point( alpha=0.7, size=1.5, shape=21, stroke=0.5) + theme(aspect.ratio = 1) + 
          scale_colour_gradient(low="gray90",high="red",guide = "colourbar", name=g1) + theme_bw()  +  
          scale_fill_gradient(low="gray90",high="blue",guide = "colourbar",name=g2) + theme(aspect.ratio = 1, legend.position="bottom") + xlab("tSNE_2") + ylab("tSNE_2"))
}

plot1geneSeuratTSNE = function(x,g1){
  print(ggplot(data.frame(x@dr$tsne@cell.embeddings), aes(x=tSNE_1,y=tSNE_2, color=as.numeric(x@data[g1,]))) +
          geom_point( alpha=0.7, size=1.5, shape=21, stroke=0.5) + theme(aspect.ratio = 1) + 
          scale_colour_gradient(low="gray90",high="red",guide = "colourbar", name=g1) + theme_bw()  + theme(aspect.ratio = 1, legend.position="bottom") + xlab("tSNE_2") + ylab("tSNE_2"))
}

plotMetaDataSeuratPCA = function(x,m){
  
  ggplot(data.frame(x@dr$pca@cell.embeddings), aes(x=PC1,y=PC2, color=x@meta.data[,m])) + geom_point( alpha=0.6, size=1) + theme(aspect.ratio = 1) + theme_bw()  +  guides(color=guide_legend(title=m)) + theme(aspect.ratio = 1, legend.position="bottom") + xlab("PC 1") + ylab("PC 2")
  
}


myFeaturePlot = function(x,m){
  
  ggplot(data.frame(x@dr$tsne@cell.embeddings), aes(x=tSNE_1,y=tSNE_2, color=x@meta.data[,m])) + geom_point( alpha=0.6, size=1)  + scale_colour_gradient(low="gray",high="blue") + theme(aspect.ratio = 1) + theme_bw()  +  guides(color=guide_legend(title=m)) + theme(aspect.ratio = 1, legend.position="bottom") + xlab("tSNE 1") + ylab("tSNE 2")
  
}

myFeaturePlotAUC = function(x,m,cells_AUC){
  
  ggplot(data.frame(x@dr$tsne@cell.embeddings), aes(x=tSNE_1,y=tSNE_2, color=assays(cells_AUC)$AUC[m,])) + geom_point( alpha=0.6, size=1)  + scale_colour_gradient(low="gray",high="blue") + theme(aspect.ratio = 1) + theme_bw()  +  guides(color=guide_legend(title=m)) + theme(aspect.ratio = 1, legend.position="bottom") + xlab("tSNE 1") + ylab("tSNE 2")
  
}


gg_color_hue <- function(n,alpha=1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#
library(dplyr)
library(tidyr)

myDotPlot = function (object, genes.plot, cols.use = c("lightgrey", "blue"), 
                      col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                      scale.by = "radius", scale.min = NA, scale.max = NA, group.by, 
                      plot.legend = FALSE, do.return = FALSE, x.lab.rot = FALSE) 
{
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  if (!missing(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
  colnames(x = data.to.plot) <- genes.plot
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- object@ident
  data.to.plot <- data.to.plot %>% gather(key = genes.plot, 
                                          value = expression, -c(cell, id))
  data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% 
    summarize(avg.exp = mean(expm1(x = expression)), pct.exp = Seurat:::PercentAbove(x = expression, 
                                                                                     threshold = 0))
  data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>% 
    dplyr::mutate(avg.exp.scale = scale(x = avg.exp,center=F,scale = T)) %>% dplyr::mutate(avg.exp.scale = MinMax(data = avg.exp.scale, 
                                                                                                    max = col.max, min = col.min))
  data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot, 
                                    levels = rev(x = genes.plot))
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, 
                                                 y = id)) + geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) + 
    scale.func(range = c(0, dot.scale), limits = c(scale.min, 
                                                   scale.max)) + theme(axis.title.x = element_blank(), 
                                                                       axis.title.y = element_blank())
  if (length(x = cols.use) == 1) {
    p <- p + scale_color_distiller(palette = cols.use)
  }
  else {
    p <- p + scale_color_gradient(low = cols.use[1], high = cols.use[2])
  }
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, 
                                              vjust = 0.5))
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}

#

