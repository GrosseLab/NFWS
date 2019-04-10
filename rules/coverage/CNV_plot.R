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

cat("remove chrY")

dt <- dt[Chromosome != "chrY",]

dt[,Feature := sub("-.*", "", dt[,FeatureRegionName])]

mean_Cov_perSample<-dt[,mean(Coverage), by=Sample]
setkey(mean_Cov_perSample, Sample)
setkey(dt, FeatureRegionName, Position, Sample)
dt[,scaledCoverage := Coverage / mean_Cov_perSample[, V1] * mean(Coverage)]

queomean<-function(x){tmp<-x[-order(abs(x-mean(x)), decreasing=T)[1]];return(prod(tmp)^(1/length(tmp)))}

median_scaledCov_perPosition<-dt[,median(scaledCoverage), by=.(FeatureRegionName, Position)]
setkey(median_scaledCov_perPosition, FeatureRegionName, Position)
setkey(dt, Sample, FeatureRegionName, Position)
dt[,relScaledCoverage := 2* scaledCoverage / median_scaledCov_perPosition[, V1] ]

coverage_hist<-dt[,list(sum=length(relScaledCoverage)) , by=round(scaledCoverage, digits = -1)+1]
N<-sum(coverage_hist[,sum])
coverage_hist[, sum := sum/N]
setkey(coverage_hist, round)
dt[,coverageWeight := coverage_hist[.(round(scaledCoverage, digits = -1)+1),sum]]

nnquantile_perCoverage<-dt[,list(lower=quantile(relScaledCoverage, 0.005 * unique(coverageWeight), na.rm=T, type=7), upper=quantile(relScaledCoverage, 1 - 0.005 * unique(coverageWeight), na.rm=T, type=7)), by=round(scaledCoverage, digits = -1)+1]
setkey(nnquantile_perCoverage, round)

median_scaledCov_perPosition[,lower := nnquantile_perCoverage[.(round(V1, digits = -1)+1), lower]]
median_scaledCov_perPosition[,upper := nnquantile_perCoverage[.(round(V1, digits = -1)+1), upper]]
setkey(median_scaledCov_perPosition, FeatureRegionName, Position)
setkey(dt, Sample, FeatureRegionName, Position)
dt[,lower := median_scaledCov_perPosition[, lower] ]
dt[,upper := median_scaledCov_perPosition[, upper] ]

split_table <- split(dt, dt[,Feature])
rm(dt)

cat("multithreaded adjustments")

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

cat("multithreaded plotting")

lapply(split_table, function(f_table){
  f <- f_table[1,Feature]

  cat(paste("Plotting", f_table[,Name][1], f, "\n"))
  features <- split(f_table, f=f_table[,Region])
  features <- lapply(features, function(l){split(l, l[,PlotRegion])})
  features <- unlist(features, recursive=F)

  nd<-function(x, top=10){return(ceiling(x/ceiling(x/top)))}

  features2 <- list()

  plots_per_page <- nd(length(features))

  for (i in seq(1,length(features), plots_per_page)) {
    features2 <- c(features2, list(do.call(rbind, features[i:min(i+plots_per_page-1, length(features))])))
  }

  plots<-lapply(features2, function(currfeature){
        p <- ggplot( data = currfeature, aes( x = Position, y = relScaledCoverage, colour = Sample))
        p <- p + geom_line()
        p <- p + geom_line(aes(y = upper, colour = c("0.95 Quantil")), colour = "grey70", size = 0.5)
        p <- p + geom_line(aes(y = lower, colour = c("0.05 Quantil")), colour = "grey70", size = 0.5)

        p <- p + expand_limits(y = c(0, 4))
        # p <- p + geom_vline( xintercept = 0.5, colour="grey", linetype = "longdash" )
        # p <- p + geom_vline( aes( xintercept = End - Start - 20 - 20 + 0.5 ), colour="grey", linetype = "longdash" )

        p <- p + facet_wrap( ~ Label, scales = 'free_x', ncol=1 )
        
        p <- p + labs(title = paste0('CNV-Analyse von ', f, '(', currfeature[1,Name], ') in ', snakemake@wildcards$ref, ', Lauf: ', snakemake@config$run), x='Position in Exon', y='relative normalisierte Coverage [erwartete Kopienanzahl]')
        p <- p + theme(plot.title = element_text(size = rel(0.8)))
        return(p)
  })

  pdf(paste0('plots/CNV/', snakemake@wildcards$mapper, snakemake@wildcards$parameterSet, '/', snakemake@wildcards$ref, '/', snakemake@wildcards$kind, '/', f_table[,Name][1], '-', f, '.pdf'), paper="a4", height=11.692, width=8.267)
    lapply(plots, print)
  dev.off()
  return(NULL)
})#, mc.cores = {threads})
