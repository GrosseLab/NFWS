#! /usr/bin/perl

use Getopt::Std;
use File::Basename;
use warnings;
use strict;

my %options=();
getopts("hr:n:m:sf:b:dg:", \%options);

eval {

	my %ref;
	($ref{"filename"}, $ref{"dir"}, $ref{"suffix"}) = fileparse($options{r}, qr/\.[^.]*/);
	
	my %reads;
	($reads{"filename"}, $reads{"dir"}, $reads{"suffix"}) = fileparse($options{n}, qr/\.[^.]*/);
	
	my $base = "./";
	if (exists $options{b}) {$base = $options{b} . "/";}
	
	our $start = 0;
	if (exists $options{f}) {$start = $options{f};}
	our $step = 0;

	my $readgroup = "RGID=foo RGLB=bar RGSM=Sample";
	if (exists $options{g}) {$readgroup = $options{g};}

	# unzip files if needed
		if ($reads{"suffix"} eq ".gz") {
			if ($start <= 1) {
				execprint("pigz -p 60 -d \"$options{n}\"");
				if (exists $options{m}) {
					execprint("pigz -p 60 -d \"$options{m}\"");
					$options{m} = substr($options{m}, 0, -3);
				}
			}
			$options{n} = substr($options{n}, 0, -3);
			($reads{"filename"}, $reads{"dir"}, $reads{"suffix"}) = fileparse($options{n}, qr/\.[^.]*/);
		}
		elsif ($reads{"suffix"} eq ".bz2") {
			if ($start <= 1) {
				execprint("bunzip2 \"$options{n}\"");
				if (exists $options{m}) {
					execprint("bunzip2 \"$options{m}\"");
					$options{m} = substr($options{m}, 0, -4);
				}
			}
			$options{n} = substr($options{n}, 0, -4);
			($reads{"filename"}, $reads{"dir"}, $reads{"suffix"}) = fileparse($options{n}, qr/\.[^.]*/);
		}

	# creating mapping and other indexes
		my $mapping_index = $ref{'dir'} . $ref{'filename'} . ".idx";
		execprint("segemehl.x -d \"$options{r}\" -x \"$mapping_index\"") if (! -e $mapping_index);
		
		execprint("samtools faidx $options{r}") if (! -e $options{r} . ".fai");

		my $ref_dict = $ref{'dir'} . $ref{'filename'} . ".dict";
		execprint("java -jar /usr/local/bin/picard-tools-1.119/CreateSequenceDictionary.jar R=$options{r} O=$ref_dict") if (! -e $ref_dict);

		execprint("bwa index -a bwtsw -p \"$ref{'dir'}$ref{'filename'}\" \"$options{r}\"") if (! -e $ref{'dir'} . $ref{'filename'} . ".sa");

	# Quality control
		system('fastqc "' . $options{n} . '" -o ' . $base . 'results/qa/fastqc/ &') if $start <= $step;
		system('fastqc "' . $options{m} . '" -o ' . $base . 'results/qa/fastqc/ &') if $start <= $step and exists $options{m};

	$step++; # 1: Mapping with segemehl
		my %sam_file;
		($sam_file{"filename"}, $sam_file{"dir"}, $sam_file{"suffix"}) = fileparse($base . "results/mapping/segemehl/segemehl_" . $reads{"filename"} . "_on_" . $ref{"filename"} . "_D1_H0.sam", qr/\.[^.]*/);
		my $unmatched_file = $sam_file{'dir'} . $sam_file{'filename'} . "_unmatched.fastq";
		my $query = "-q \"$options{n}\""; $query = $query . " -p \"$options{m}\"" if (exists $options{m});
		# execprint("segemehl.x -i \"$mapping_index\" -d \"$options{r}\" $query -s -D 1 -H 0 -t 60 -u \"$unmatched_file\" > \"$sam_file{'dir'}$sam_file{'filename'}$sam_file{'suffix'}\"");
		#1: Mapping with BWA
		($sam_file{"filename"}, $sam_file{"dir"}, $sam_file{"suffix"}) = fileparse($base . "results/mapping/bwa/bwa_" . $reads{"filename"} . "_on_" . $ref{"filename"} . "_M.sam", qr/\.[^.]*/);
		$query = "\"$options{n}\""; $query = $query . " \"$options{m}\"" if (exists $options{m});
		execprint("bwa mem -t 60 -M \"$ref{'dir'}$ref{'filename'}\" $query > \"$sam_file{'dir'}$sam_file{'filename'}$sam_file{'suffix'}\"");

		system('pigz -p 20 "' . $options{n} . '" &') if $start <= $step;
		# system('pigz -p 20 "' . $unmatched_file . '" &') if $start <= $step;
		system('pigz -p 20 "' . $options{m} . '" &') if $start <= $step and exists $options{m};

	$step++; # 2: convert sam to bam
		my $bam_file = $sam_file{"dir"} . $sam_file{"filename"} . ".bam";
		execprint("samtools view -@ 60 -b -o \"$bam_file\" \"$sam_file{'dir'}$sam_file{'filename'}$sam_file{'suffix'}\"");
		execprint("rm \"$sam_file{'dir'}$sam_file{'filename'}$sam_file{'suffix'}\"") if (! exists $options{d});

	$step++; # 3: sort bam
		my $sorted_bam_file = $sam_file{"dir"} . $sam_file{"filename"} . "_sort.bam";
		execprint("samtools sort -m 512M -@ 60 -o \"$sorted_bam_file\" -T tmp_samtools_sort \"$bam_file\"");
		execprint("rm \"$bam_file\"") if (! exists $options{d});

	$step++; # 4: add MappingQuality -> Fabian im November

	$step++; # 5: mark duplicates with Picard
		my $dedup_sorted_bam_file = $sam_file{"dir"} . $sam_file{"filename"} . "_sort_dedup.bam";
		execprint("java -Xmx40g -jar /usr/local/bin/picard-tools-1.119/MarkDuplicates.jar INPUT=$sorted_bam_file OUTPUT=$dedup_sorted_bam_file METRICS_FILE=$base/results/qa/picard_duplication_report.txt READ_NAME_REGEX=\"HWI-[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\" MAX_RECORDS_IN_RAM=5000000");
		execprint("rm \"$sorted_bam_file\"") if (! exists $options{d});
	
	$step++; # 6: add ReadGroups
		my $RG_dedup_sorted_bam_file = $sam_file{"dir"} . $sam_file{"filename"} . "_sort_dedup_RG.bam";
		execprint("java -jar /usr/local/bin/picard-tools-1.119/AddOrReplaceReadGroups.jar I=$dedup_sorted_bam_file O=$RG_dedup_sorted_bam_file $readgroup RGPL=illumina RGPU=na");
		execprint("rm \"$dedup_sorted_bam_file\"") if (! exists $options{d});

	$step++; # 7: index bam
		execprint("samtools index \"$RG_dedup_sorted_bam_file\"");

	$step++; # 8: Indel Realignment with GATK
		my $realigned_RG_dedup_sorted_bam_file = $sam_file{"dir"} . $sam_file{"filename"} . "_sort_dedup_RG_realign.bam";
		my $intervalList = substr(`mktemp intervalList.XXXXX.intervals`, 0, -1);
		execprint("java -Xmx2g -jar /usr/local/bin/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R \"$options{r}\" -I \"$RG_dedup_sorted_bam_file\" -o \"$intervalList\" -nt 60");
		execprint("java -Xmx4g -jar /usr/local/bin/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -R \"$options{r}\" -I \"$RG_dedup_sorted_bam_file\" -targetIntervals \"$intervalList\" -o \"$realigned_RG_dedup_sorted_bam_file\"");
		execprint("rm \"$intervalList\" \"$RG_dedup_sorted_bam_file\"*") if (! exists $options{d});
	
	$step++; # 9: Base recalibration

	$step++; # 10: statistics
		my $stat_file = $realigned_RG_dedup_sorted_bam_file . "_stat.txt";
		system('samtools idxstats "' . $realigned_RG_dedup_sorted_bam_file . '" > "' . $stat_file . '"') if $start <= $step;;

	$step++; # 11: RR compression

	$step++; # 12: call SNPs
		my $gvcf_file = $base . "results/snp/gatk/" . $sam_file{"filename"} . "_variants.g.vcf";
		execprint("java -jar /usr/local/bin/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T HaplotypeCaller -R $options{r} -I $realigned_RG_dedup_sorted_bam_file --emitRefConfidence GVCF --variant_index_type LINEAR -nct 60 --variant_index_parameter 128000 -o $gvcf_file");

	$step++; # 13: convert g.vcf file to single sample vcf file
		my $vcf_file = $base . "results/snp/gatk/" . $sam_file{"filename"} . "_variants.vcf";
		execprint("java -jar /usr/local/bin/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $options{r} --variant $gvcf_file -o $vcf_file");

	$step++; # 14: analyse SNPs with Analysis-Pipeline
		execprint("perl ~/NGS-scripts/DefaultPipelines/AnalyseSNPs.pl -i $vcf_file -7");
		execprint("rm $vcf_file");
		execprint("rm " . $vcf_file . ".idx");

	system("echo 'SNP-Calling of $options{n} on $options{r} is done'|mail -s 'SNP-Calling is done' martin.porsch\@informatik.uni-halle.de") if (! exists $options{s});
};
system("echo 'SNP-Calling of $options{n} on $options{r} failed with:\n$@' | mail -s 'SNP-Calling error' martin.porsch\@informatik.uni-halle.de") if ($@);

sub execprint {
	print STDERR "$main::step: " . $_[0] . "\n";
	print `$_[0]` if ($main::start <= $main::step);
	die $? if $?;
}