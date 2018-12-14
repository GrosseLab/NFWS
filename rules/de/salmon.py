import tempfile
from os import system

if len(snakemake.input.bams) == 1:
	bam_file = snakemake.input.bams[0]
else:
	bam_file = tempfile.mkstemp(prefix="./", suffix=".bam")[1]
	system("samtools merge -f " + bam_file + " " + snakemake.input.bams)

system("salmon quant --threads " + str(snakemake.threads) + " --targets " + snakemake.input.transcriptome + " --alignments " + bam_file + " --geneMap " + snakemake.input.gtf + " " + " ".join(snakemake.params))

if len(snakemake.input.bams) > 1:
	os.remove(bam_file)
