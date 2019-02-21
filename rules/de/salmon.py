import tempfile
from os import system
from os import remove

if len(snakemake.input.bams) == 1:
	bam_file = snakemake.input.bams[0]
else:
	bam_file = tempfile.mkstemp(prefix="./", suffix=".bam")[1]
	system("samtools merge -f " + bam_file + " " + " ".join(snakemake.input.bams))

system("salmon quant --threads " + str(snakemake.threads) + " --targets " + snakemake.input.transcriptome + " --alignments " + bam_file + " " + " ".join(snakemake.params))

if len(snakemake.input.bams) > 1:
	remove(bam_file)
