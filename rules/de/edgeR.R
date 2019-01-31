library(edgeR)
# library(affy)
# install.packages("pheatmap")
library(pheatmap)
library(ggplot2)
# library(tximport)

saveRDS(snakemake, file="debug_edger.rds")

plotPCA <- function(exp_mat, groups, do.legend=T, log=T, do.MDS=F,do.ggplot=T,epsilon=1, plot_label = F, plot_chull=TRUE,plot_geom_path=FALSE,plot_title=""){
  
  if(!is.factor(groups)){
    groups <- factor(groups)
  }
  
  group_levels <- droplevels(groups)
  rep_indices <- split(seq_along(groups), f = group_levels)
  
  nCols <- length(unique(group_levels))
  
  labelColors <- suppressWarnings(RColorBrewer::brewer.pal(name = "Set1",n = nCols))
  labelColors <- scales::alpha(colour = labelColors, alpha = 0.8)
  
  j <- 0
  rep_groups <- lapply(rep_indices, function(i){
    j <<- j +1
    cbind(colnames(exp_mat)[i],names(rep_indices)[j], labelColors[j])
  })
  
  sample_group <- do.call("rbind",rep_groups)
  colnames(sample_group) <- c("Sample","Group","Color")
  rownames(sample_group) <- sample_group[,1]
  
  if(log){
    exp_mat <- log10(exp_mat + epsilon)
  }
  
  x.lab <- c()
  y.lab <- c()
  
  x <- c()
  y <- c()
  if(do.MDS){
    tmp.dist <- dist(t(exp_mat))
    fit <- cmdscale(tmp.dist,eig=TRUE, k=2) # k is the number of dim
    pca <- list()
    x <- fit$points[,1]
    y <- fit$points[,2]
    
    x.lab <- "Coordinate 1"
    y.lab <- "Coordinate 2"
  }else{
    pca <- prcomp(t(exp_mat))
    x <- pca$x[,1] 
    y <- pca$x[,2]
    
    s <- summary(pca)
    x.lab <- paste0("PC1 (",round(s[[6]][2,1]*100,1),"%)")
    y.lab <-paste0("PC2 (",round(s[[6]][2,2]*100,1),"%)")
  }
  
  res.list <- list()
  
  plot.data <- data.frame(x=x, y=y,Groups = sample_group[names(x),2], Colors = sample_group[names(x),3])
  #plot.data$Groups <- factor(plot.data$Groups, labels = plot.data$Groups)
  
  if(do.ggplot){
    
    # calculate minimal convex hull
    if(plot_chull) chulls <- plyr::ddply(plot.data, plyr::.(Groups), function(df) df[chull(df$x, df$y), ])
    
    pcaPlot <- ggplot2::ggplot(plot.data, ggplot2::aes(x=x, y=y, group=Groups, col=Groups)) + ggplot2::geom_point(size=2)  + ggplot2::labs( x=x.lab,  y=y.lab)  

    if(plot_geom_path){
      plot_chull <- FALSE
      pcaPlot <- pcaPlot + geom_path(arrow = arrow(angle = 30, length = unit(0.2, "inches"), ends = "last", type='closed')) # ,type = "open"))
    }

    if(plot_chull){
      pcaPlot <- pcaPlot + ggplot2::geom_polygon(data = chulls, ggplot2::aes(x=x, y=y, fill=Groups), alpha = 0.2) 
    }
    
    if(plot_label){
      pcaPlot <- pcaPlot + ggplot2::geom_text(ggplot2::aes(label=rownames(plot.data)),hjust=0, vjust=0)
    }
    
    if(plot_title!=""){
      pcaPlot <- pcaPlot + labs(title=plot_title) + theme(plot.title = element_text(hjust = 0.5))
    }


    res.list <- list(plot=pcaPlot, plot.data = plot.data )
  }else{
    plot(plot.data$x, plot.data$y, col = plot.data$Colors, xlab = x.lab, ylab = y.lab, pch=19)
    if(plot_label){
      text(x = plot.data$x, y = plot.data$y, labels = rownames(plot.data))
    }
    if(do.legend){
      legend("topleft",legend = unique(plot.data$Groups),col = unique(plot.data$Colors),bty = "n",pch = 19,ncol=2,pt.cex = 2)    
    }
    res.list <- list(plot.data = plot.data)
  }
  
  return(res.list)
}

counts <- read.table(snakemake@input[[1]])
length <- read.table(snakemake@input[[2]])

rep_list <- snakemake@config[["replicates"]]
rep_names <- names(rep_list)
rep_dict <- rep(rep_names, lapply(rep_list, length))
names(rep_dict) <- gsub("^(?=[0-9])", "X", gsub("-", ".", do.call(c, rep_list[rep_names])), perl=T)
group <- factor(rep_dict[colnames(counts)], levels = rep_names)

normMat <- length
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(counts/normMat)) + log(colSums(counts/normMat))
y <- DGEList(counts = counts, group	= group)
y$offset <- t(t(log(normMat)) + o)

y <- calcNormFactors(y, logratioTrim=0.499, sumTrim=0.0, doWeighting=T)

rpkm <- rpkm(y, gene.length=apply(length, 1, mean))
tpm <- scale(rpkm, center=F, scale=colSums(rpkm))*1e6

# mean over replicates
summarize_replicates <- function(exp_mat, groups, method=mean, changeColNames = T){
	if(!is.factor(groups)){
		groups <- factor(groups)
	}
	
	if(length(groups) != ncol(exp_mat)){
		stop("Length of your group vector is not equal to the number of columns in your expression matrix.")
	}
	
	group_levels <- droplevels(groups)
	rep_indices <- split(seq_along(groups), f = group_levels)
	sum_exp_mat <- sapply(rep_indices, function(i){
		apply(exp_mat[,i,drop=FALSE], 1, method)
	})
	if(changeColNames){
		colnames(sum_exp_mat) <- levels(groups)
	}else{
		new_name_idc <- sapply(rep_indices, `[`, 1)
		colnames(sum_exp_mat) <- colnames(exp_mat[,new_name_idc])
	}
	
	return(sum_exp_mat)
}
tpm_avrg<-summarize_replicates(tpm, group)
rpkm_avrg<-summarize_replicates(rpkm, group)

# write out tables
write.csv(tpm_avrg, file=snakemake@output[["tpms"]], quote=F)
write.csv(rpkm_avrg, file=snakemake@output[["rpkms"]], quote=F)

# write out plots
pdf(snakemake@output[["plot1"]])
	plot(hclust(dist(t(tpm))), xlab="Sample")
dev.off()

pdf(snakemake@output[["plot2"]])
	plot(hclust(as.dist(1-cor(tpm, method="spearman")), method="average"), xlab="Sample")
dev.off()

p <- plotPCA(tpm, group, do.MDS=T, do.ggplot=T, log = T, plot_label = T)$plot#, labels=unique(group)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(snakemake@output[["plot3"]], p, width=6, height=6)

pdf(snakemake@output[["plot4"]])
	plotDensities(log(tpm+1), log=T, legend="topright")
dev.off()

pdf(snakemake@output[["plot5"]])
	boxplot(log(tpm+1))
dev.off()

pdf(snakemake@output[["plot6"]])
	boxplot(t(t(log10(counts+1))-apply(log10(counts+1), 1, median)), outline=F)
dev.off()

pheatmap(cor(tpm), file=snakemake@output[["plot7"]])
pheatmap(cor(log(tpm+1)), file=snakemake@output[["plot8"]])
pheatmap(cor(tpm, method = "spearman"), file=snakemake@output[["plot9"]])

pdf(snakemake@output[["plot_libSize"]])
	barplot(y$samples$lib.size)
dev.off()
