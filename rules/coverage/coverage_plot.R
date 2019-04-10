library(ggplot2)
library(data.table)
library(parallel)

files <- unlist(snakemake@input)
labels <- sub( '_20_20.cov', "", basename(files)) # extract sample names
coverages<-list()
for (i in seq_along(files)) {
  coverages[[i]] <- fread(files[i], sep="\t")
  coverages[[i]] <- cbind( labels[i], coverages[[i]])
  setnames(coverages[[i]], c("Sample", "Chromosome", "Start", "End", "FeatureRegionName", "Position", "Coverage"))
}

dt <- do.call('rbind', coverages)
rm(coverages)

dt[,Feature := sub("-.*", "", dt[,FeatureRegionName])]

split_table <- split(dt, dt[,Feature])
rm(dt)

cat("multithreaded updates")

split_table <- lapply(split_table, function(f_table){
  f_table[,RegionName := sub(".*?-", "", f_table[,FeatureRegionName])]
  f_table[,Region := as.numeric(sub("-.*", "", f_table[,RegionName]))]
  f_table[,Name := sub(".*-", "", f_table[,FeatureRegionName])]
  f_table[,PlotRegion := (Position %/% 500) + 1]

  f_table[,Position := f_table[,Position - 20]]
  f_table[,LocationText := paste0(f_table[,Chromosome], ":", format(f_table[,Start] + 20 - 1, big.mark = ".", decimal.mark=","), "-", format(f_table[,End] - 20, big.mark = ".", decimal.mark=","))]
  # f_table[,Exon := "Exon"]
  # f_table[,Abschnitt := "Abschnitt"]
  f_table[,Label := paste0("Exon ", Region, " (", LocationText, "), Abschnitt ", PlotRegion)]
  setkey(f_table, Region, Position, Sample)
  return(f_table)
})#, mc.cores = {threads})

cat("starting plotting")

lapply(split_table, function(f_table){

  f <- f_table[1,Feature]

  cat(paste("Plotting", f_table[,Name][1], f))
  features <- split(f_table, f_table[,Region])
  features <- lapply(features, function(l){split(l, l[,PlotRegion])})
  features <- unlist(features, recursive=F)

  nd<-function(x, top=10){return(ceiling(x/ceiling(x/top)))}

  features2 <- list()

  plots_per_page <- nd(length(features))

  for (i in seq(1,length(features), plots_per_page)) {
    features2 <- c(features2, list(do.call(rbind, features[i:min(i+plots_per_page-1, length(features))])))
  }

  plots<-lapply(features2, function(features){
        p <- ggplot( data = features, aes( x = Position, y = Coverage, colour = Sample ))
        p <- p + geom_line() + geom_hline( aes(yintercept = 100)) 

        # p <- p + geom_vline( xintercept = 0.5, colour="grey", linetype = "longdash" )
        # p <- p + geom_vline( aes( xintercept = End - Start - {left} - {right} + 0.5 ), colour="grey", linetype = "longdash" )
        # p <- p + geom_vline( aes( xintercept = End - Start - 20 - 20 + 0.5 ), colour="grey", linetype = "longdash" )

        p <- p + facet_wrap( ~ Label, scales = 'free_x', ncol=1 )
        
        p <- p + labs(title = paste0('Abdeckung von ', f, '(', features[1,Name], ') in ', snakemake@wildcards$ref, ', Lauf: ', snakemake@config$run), x='Position in Exon' )
        p <- p + theme(plot.title = element_text(size = rel(0.8)))
        return(p)
  })

  pdf(paste0('plots/coverage/' , snakemake@wildcards$mapper, snakemake@wildcards$parameterSet, '/', snakemake@wildcards$ref, '/', snakemake@wildcards$kind, '/', f_table[,Name][1], '-', f, '.pdf'), paper="a4", height=11.692, width=8.267)
    lapply(plots, print)
    lapply(plots, function(p){
      print(p + scale_y_log10(limits = c(1, 20000), breaks = 10^(1:4))
    )}) 
  dev.off()
  return(NULL)
})#, mc.cores = {threads})
