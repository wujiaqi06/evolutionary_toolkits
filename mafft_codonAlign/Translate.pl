#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;

my %opts;
getopts('i:o:t:', \%opts);
my $input_folder = $opts{'i'} or die "use: perl $0 -i Genome_folder -o output_folder -t codon_table\n";
my $output_folder = $opts{'o'} or die "use: perl $0 -i Genome_folder -o output_folder -t codon_table\n";
my $codon_table = $opts{'t'} or die "use: perl $0 -i Genome_folder -o output_folder -t codon_table\n";


open TABLE, $codon_table or die "Cannot open Standard.txt: $!\n";
my %table;
while (<TABLE>){
	chomp;
	my @lines = split;
	$table{$lines[0]} = $lines[1];
}

chdir "./$input_folder" or die "Cannot enter cds folder: $!\n";
my @cds = glob ("*.fas");
chdir "../";
#print "@cds\n";
mkdir "./$output_folder" or die "Cannot make folder $output_folder: $!\n";

foreach my $seq (@cds){
	my %cds = &read_fasta("./$input_folder/$seq");
	my $seq0;
	if ($seq =~ /(.*)\.fas/){
		$seq0 = $1;
	}
	open PEP, ">./$output_folder/$seq0.pep.fas" or die "Can not make file $seq.pep.fas\n";
	my @cds_name = sort keys %cds;
	foreach my $i (@cds_name){
		my @cds_seq = split(//,$cds{$i});
		#my $codon_pos = $cds_seq[3].$cds_seq[4].$cds_seq[5];
		#print "$codon_pos\n";
		my $pep_length = int (($#cds_seq + 1)/3);
		my $pep_seq;
		foreach my $j (0..($pep_length-1)){
			my $codon_pos = $cds_seq[3*$j].$cds_seq[3*$j+1].$cds_seq[3*$j+2];
			$codon_pos = uc $codon_pos;
			if (exists $table{$codon_pos}){
				$pep_seq .= $table{$codon_pos};
			} elsif ($codon_pos eq "---") {
				#print "$codon_pos\n";
				$pep_seq .= "-";
			}	else {
				$pep_seq .= "-"
			}
		}
		my @pep_site = split (//, $pep_seq);
		if ($pep_site[-1] eq "*"){
			pop (@pep_site);
		}
		print PEP ">$i\n";
		print PEP @pep_site;
		print PEP "\n";
	}
}




sub read_fasta{
	open READ_SEQ, "@_" or die "Can not read Sequence file: $!\n";
	open SEQ_INFO, ">>../Sequence_Info.txt" or die "Can not read Sequence file: $!\n";
	my %sequences;
	my $id;
	my $count = 0;
	my $seq_length;
	while (<READ_SEQ>){
	chomp;
	s/\r//;
	if ($_ =~ /^>(?<seq_name>.+$)/){
		$id  = $+{seq_name};
		$count += 1;
	}else{
		s/-//g;
		$sequences{$id} .= $_;
		$seq_length = length ($sequences{$id});
		}
	}
	print "Successfuly read @_, in total $count sequences with $seq_length bp.\n";
	print SEQ_INFO "Successfuly read @_, in total $count sequences with $seq_length bp.\n";
	close SEQ_INFO;
	#print "\n";
	%sequences;
}


