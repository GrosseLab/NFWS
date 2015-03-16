#! /usr/bin/perl

use Getopt::Std;
use File::Basename;
use warnings;
use strict;

my %options=();
getopts("hi:g:78s", \%options);

eval {
	my %vcf;
	($vcf{"filename"}, $vcf{"dir"}, $vcf{"suffix"}) = fileparse($options{i}, qr/\.[^.]*/);

	our $start = 0;
	if (exists $options{g}) {$start = $options{g};}
	our $step = 0;

	
	my $snp_db;
	my $ref_id;
	if (exists $options{7}) {
		$snp_db = "/data/reference/dbSNP_hg19_00-All.vcf";
		$ref_id = "hg19";
	}
	elsif (exists $options{8}) {
		 $snp_db = "/data/reference/dbSNP_hg38_00-All.vcf";
		 $ref_id = "hg38";
	}
	else {
		die "specify reference version, i.e. -7 for hg19, -8 for hg38";
	}

	$step++; # 1:  annotate SNPs with dbSNP-IDs
		my $ann_vcf_file = $vcf{"dir"} . $vcf{"filename"} . "_annotated.vcf";
		execprint("java -jar /usr/local/bin/snpEff-4.0/SnpSift.jar annotate $snp_db \"$vcf{dir}$vcf{filename}$vcf{suffix}\" > $ann_vcf_file");

	$step++; # 2: predict effects with SnpEff
		my $snpEff_file = $vcf{"dir"} . $vcf{"filename"} . "_annotated_snpEff.vcf";
		execprint("java -Xmx4G -jar /usr/local/bin/snpEff-4.0/snpEff.jar -c /usr/local/bin/snpEff-4.0/snpEff.config -v $ref_id $ann_vcf_file > $snpEff_file");
		execprint("rm $ann_vcf_file");

		system("echo 'SNP-Analysis of $options{i} is done'|mail -s 'SNP-Analysis is done' martin.porsch\@informatik.uni-halle.de") if (! exists $options{s});
};
system("echo 'SNP-Analysis of $options{i} failed with:\n$@' | mail -s 'SNP-Analysis error' martin.porsch\@informatik.uni-halle.de") if ($@);

sub execprint {
	print STDERR "$main::step: " . $_[0] . "\n";
	print `$_[0]` if ($main::start <= $main::step);
	die $? if $?;
}
