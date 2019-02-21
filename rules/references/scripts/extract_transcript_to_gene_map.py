import gffutils

gff = gffutils.FeatureDB(snakemake.input.annotation_db)
output_file = open(snakemake.output[0], 'w')

for gene in gff.features_of_type("gene", order_by="start"):
	for transcript in gff.children(gene, order_by='start', level=1):
		if transcript.featuretype != "exon":
			output_file.write(transcript['ID'][0] + "\t" + gene['ID'][0] + "\n")
		else:
			output_file.write(gene['ID'][0] + "\t" + gene['ID'][0] + "\n")
		
output_file.close()
