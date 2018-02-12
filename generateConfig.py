#! python3
import argparse
import csv
import datetime
import re
from os import listdir

# set up arguments
parser = argparse.ArgumentParser(description='Generate a config.jason file for a project from an Illumina *SeqOutput folder.')
parser.add_argument("path", type=str, help='path to the *SeqOutput folder')
parser.add_argument("-o", metavar="PATH", type=str, help='path to the config file')
parser.add_argument("--rules", type=str, help='path to the NFWS rules directory', default="/path/to/rules/dir")
parser.add_argument("--email", metavar="MAILS", type=str, help='your email adress', default="")
parser.add_argument("--estart", metavar="DATE", type=str, help='experiment start date', default="")
parser.add_argument("--keep", metavar="DATE", type=str, help='keep data at least until DATE', default="")
parser.add_argument("--pi", metavar="NAME", type=str, help='NAME of the PI respsonsible for the experiment', default="")
parser.add_argument("--researcher", metavar="NAME", type=str, help='NAME of the researcher respsonsible for the experiment', default="")
parser.add_argument("--institution", metavar="NAME", type=str, help='NAME of the institution respsonsible for the experiment', default="")
parser.add_argument("--text", metavar="TEXT", type=str, help='some descriptive TEXT', default="")

args = parser.parse_args()
# print(args)

# get Experiment name and sequencing date
run = ""
seq_date = ""
units_array = []
with open(args.path + "/SampleSheet.csv") as csvfile:
	samplesheet = csv.reader(csvfile)
	read_samples = False
	read_reads = False
	for row in samplesheet:
		if (len(row) > 0):
			if (row[0] == "Date"): seq_date = row[1]
			if (row[0] == "Experiment Name"): run = row[1]

# get current date
today = datetime.datetime.today().strftime("%d.%m.%Y")

# generate sample and unit dicts
units_array_tmp = listdir(args.path + "/Data/Intensities/BaseCalls")
units_array_tmp.sort(key=str.lower)
units_array = []
for unit in units_array_tmp:
	if (not "Undetermined" in unit and "fastq.gz" in unit):
		units_array.append(unit)

sample_array = list(set(map(lambda s: s[:-21], units_array)))
sample_array.sort(key=str.lower)

# generate samples and units strings
samples = "";
for i in range(0,len(sample_array)):
	if (i > 0):
		samples = samples + ",\n		"
	samples = samples + "\"" + sample_array[i] + "\": [ \"" + sample_array[i] + "\" ]"

units = ""
for i in range(0,len(sample_array)):
	if (i > 0):
		units = units + ",\n		"
	units = units + "\"" + sample_array[i] + "\": [\"" + "\", \"".join(filter(lambda x:sample_array[i] in x, units_array)) + "\" ]"

# print(run)
# print(today)
# print(units_array)
# print(sample_array)
# print(samples)
# print(units)

# write config
with open(args.o, "w") as config:
	config.write('{\n'
	'	"run": "' + run + '",\n'
	'	"machine": "MiSeq",\n'
	'	"platform": "illumina",\n'
	'	"project_management_version": "2.0",\n'
	'	"path_rules": "' + args.rules + '",\n\n'

	'	"description": {\n'
	'		"experiment_start": "' + args.estart + '",\n'
	'		"experiment_finish": "",\n'
	'		"sequencing date": "' + seq_date + '",\n'
	'		"analysis_start": "' + today + '",\n'
	'		"analysis_finish": "",\n'
	'		"keep_until": "' + args.keep + '",\n'
	'		"pi": "' + args.pi + '",\n'
	'		"researcher": "' + args.researcher + '",\n'
	'		"institution": "' + args.institution + '",\n'
	'		"emails": [ ] ,\n'
	'		"text": "' + args.text + '",\n'
	'		"publication": [""]\n'
	'	},\n\n'

	'	"ref_dir": "/data/references/human",\n'
	'	"references": {\n'
	'		"hg38": {\n'
	'			"genome": "/hg38/genome/hg38.fa",\n'
	'			"genome_annotation": "/hg38/genes/all.gtf",\n'
	'			"transcriptome": "/hg38/genes/all.gtf.transcript.fa"\n'
	'		},\n'
	'		"hg19": {\n'
	'			"genome": "/hg19/genome/hg19_ucsc_rm.fa",\n'
	'			"genome_annotation": "/hg19/genes/all.gtf",\n'
	'			"transcriptome": "/hg19/genes/all.gtf.transcript.fa"\n'
	'		}\n'
	'	},\n'
	'	"adapter_sequences": "/data/references/illumina/adapters/TruSeq3-PE.fa",\n\n'

	'	"samples": {\n'
	'		' + samples + '\n'
	'	},\n'
	'	"units": {\n'
	'		' + units + '\n'
	'	},\n\n'

	'	"parameters": {\n'
	'		"parameterSet": {\n'
	'			"program1": ["--arg", "value"]\n',
	'			"program2": ["--arg", "value"]\n'
	'		}\n'
	'	}\n'
	'}')
