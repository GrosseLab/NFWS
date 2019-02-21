library(tximport)

samples <- names(snakemake@config[["samples"]])
files <- paste0("results/de/salmon", snakemake@wildcards[["parameterSet"]], "/", snakemake@wildcards[["ref"]], "/", snakemake@wildcards[["kind"]], "/", samples, "/quant.sf")
names(files) <- samples

txi <- tximport(files, type = "salmon", txOut = TRUE, importer = read.delim)

write.table(txi$counts, file=snakemake@output[[1]], sep="\t", quote=F)
write.table(txi$length, file=snakemake@output[[2]], sep="\t", quote=F)
