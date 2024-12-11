#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;

my %opts;
getopts('i:o:', \%opts);
my $input_folder = $opts{'i'} or die "use: $0 -i Genome_folder -o output_folder\n";
my $output_folder = $opts{'o'} or die "use: $0 -i Genome_folder -o output_folder\n";

mkdir "$output_folder" or die "Can not create folder $output_folder:$!\n";
chdir "./$input_folder";
my @all_seq_file = glob '*.fas*';
foreach my $seqfile (@all_seq_file){
	my $seq_name;
	if ($seqfile =~ /(.*)\.fas(ta)?/){
		$seq_name = $1;
	}
	system ("mafft --auto $seqfile > .././$output_folder/$seq_name.align.fas");
}