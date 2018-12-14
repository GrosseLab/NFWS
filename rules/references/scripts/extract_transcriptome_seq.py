import gffutils
from Bio import SeqIO
fasta = SeqIO.to_dict(SeqIO.parse(snakemake.input.genome, "fasta"))
gff = gffutils.FeatureDB(snakemake.input.annotation_db)
output_file = open(snakemake.output[0], 'w')

for feature in gff.features_of_type(snakemake.wildcards.feature, order_by="start"):
	output_file.write(">" + feature['ID'][0] + "\n")
	for exon in gff.children(feature, featuretype='exon', order_by='start'):
		output_file.write(str(fasta[feature.seqid].seq)[exon.start-1:exon.stop])
	output_file.write("\n")

output_file.close()
