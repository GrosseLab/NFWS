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

# no Y-Chromosome since it messes up the analysis
dt<-dt[Chromosome != "chrY" & Chromosome != "chrX",]

dt[,Feature := sub("-.*", "", dt[,FeatureRegionName])]

mean_Cov_perSample<-dt[,mean(Coverage), by=Sample]
setkey(mean_Cov_perSample, Sample)
setkey(dt, FeatureRegionName, Position, Sample)
dt[,scaledCoverage := Coverage / mean_Cov_perSample[, V1] * mean(Coverage)]

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

char_sign<-function(x, l, u){ifelse(is.na(x), "?", ifelse(x>u, "+", ifelse(x<l, "-", "=")))}
dt[,Sign := char_sign(relScaledCoverage, lower, upper) ]

dt_cnv <- dt[relScaledCoverage > upper | relScaledCoverage < lower,]
dt_cnv <- split(dt_cnv[,.(Chromosome, Position, FeatureRegionName, Sign)], dt_cnv[, Sample])

sapply(names(dt_cnv), function(l_index){
	l<-dt_cnv[[l_index]]
	out<-c()
	start<-1
	i<-2
	while(i < dim(l)[1]) {
		while (i < dim(l)[1] & i-start == l[i,Position]-l[start,Position] & l[i,FeatureRegionName]==l[start,FeatureRegionName] & l[i,Sign]==l[start,Sign] & l[i,Chromosome]==l[start,Chromosome]) {
			i<-i+1
		}
		out<-rbind(out, c(l[start, .(Chromosome, FeatureRegionName, Position)], l[i-1, Position], l[start, Sign] ))
		start <- i
		i<-i+1
	}
	# dir.create(paste0("results/CNV/", snakemake@wildcards[["mapper"]], snakemake@wildcards$parameterSet, "/", snakemake@wildcards[["ref"]], "/", snakemake@wildcards[["kind"]], "/"))
	write.table(out, file = paste0("results/CNV/", snakemake@wildcards[["mapper"]], snakemake@wildcards$parameterSet, "/", snakemake@wildcards[["ref"]], "/", snakemake@wildcards[["kind"]], "/", l_index, "_CNV.txt"), row.names=F, quote=F, col.names=F, sep="\t")
})
