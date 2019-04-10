library(ggplot2)

files <- unlist(snakemake@input)
labels <- sub( '_diagnostics_20_20.hist', "", basename(files)) # extract sample names
cov <- list()

for (i in seq_along(files)) {
	cov[[i]] <- read.table(files[i], sep="\t")
	cov[[i]] <- cov[[i]][,c(-1, -4)] # drop first column
	cov[[i]] <- cbind( cov[[i]], 1 - cumsum( cov[[i]][,3] ) ) # add cummulative frequency of depth
	cov[[i]] <- cbind( labels[i], cov[[i]]) # add label column
}

df <- do.call('rbind', cov)
colnames(df) <- c('Sample', 'Depth', 'Count', 'Frequency', 'CumFrequency')

p <- ggplot( data = df, aes( x = Depth, y = CumFrequency, colour = Sample ) )
p <- p + scale_x_continuous(breaks = c(0, 30, 50, 70, 100), minor_breaks = c(150, 200), limits = c(0, 200))
p <- p + geom_hline(linetype = "5F", yintercept = unlist(lapply(split(df, df[,"Sample"]), function(df.sample){tmp<-df.sample[df.sample[,"Depth"]>=100,]; max(tmp[,"CumFrequency"])})))
p <- p + geom_line() + ylim(0.98, 1.0)
# p <- p + geom_vline( xintercept = seq(30,100,10), linetype = "longdash" )
p <- p + ggtitle(paste0("Coverage distribution in run ", snakemake@config$run, " using ", snakemake@wildcards$ref)) + xlab("Depth of coverage") 
p <- p + ylab( sprintf("Fraction of target bases with a coverage >= depth" ) )
ggsave( file = snakemake@output[[1]], plot = p, width = 8.267, height = 11.692, title = "Coverage distribution of target regions" )
