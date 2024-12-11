#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
# Usage example: perl RemoveGap.pl -i 43mammal_pro_align -o DelAmbSites_43mammal_pro_align -d 0.3 -t phy
# Usage example2: perl RemoveGap.pl -i 43mammal_pro_align -o DelAmbSites_43mammal_pro_align_phy -d 0.3 -t fas
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;

my %opts;
getopts('i:o:r:', \%opts);
my $input_folder = $opts{'i'} or die "use: perl $0 -i Genome_folder -o output_folder -r reference_seq_name";
my $output_folder = $opts{'o'} or die "use: perl $0 -i Genome_folder -o output_folder -r reference_seq_name";
my $reference_seq_name = $opts{'r'} or die "use: perl $0 -i Genome_folder -o output_folder -r reference_seq_name";

mkdir "$output_folder" or die "Can not create folder $output_folder:$!\n";
chdir "./$input_folder";
my @all_seq_file = glob '*.fas';
chdir "../";

open TABLE, "Standard.txt" or die "Cannot open Standard.txt: $!\n";
my %table;
while (<TABLE>){
	chomp;
	my @lines = split;
	$table{$lines[0]} = $lines[1];
}

foreach my $seq (@all_seq_file){
	my %read_sequences = &read_fasta("./$input_folder/$seq");
	my $ref_seq = $read_sequences{$reference_seq_name};
	my @ref_seq = split (//, $ref_seq);
	my @ref_seq_loc;
	foreach my $i (0..$#ref_seq){
		if ($ref_seq[$i] ne "-"){
			push @ref_seq_loc, $i;
		}
	}
	open OUT, ">./$output_folder/$seq" or die $!;
	foreach my $i (sort keys %read_sequences){
		my $seq = $read_sequences{$i};
		my @seq = split (//, $seq);
		my @new_seq = @seq[@ref_seq_loc];
		#print "@new_seq\n";
		my $new_seq = join ("", @new_seq);
		print OUT ">$i\n$new_seq\n";
	}
	close OUT;
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

sub translate{
	my $input = $_[0];
	#print "$input\n";
	my @cds_seq = split (//,$input);
	#print "cds length". @cds_seq . "\n";
	my $pep_seq;
	my $pep_length = int (($#cds_seq+1)/3);
	foreach my $j (0..($pep_length-1)){
		my $codon_pos = $cds_seq[3*$j].$cds_seq[3*$j+1].$cds_seq[3*$j+2];
		if (exists $table{$codon_pos}){
			$pep_seq .= $table{$codon_pos};
		} else {
			$pep_seq .= "-";
			#print "Codon is $codon_pos\n";
		}
	}
	my @pep_site = split (//, $pep_seq);
	# if ($pep_site[-1] eq "*"){
	# 	pop (@pep_site);
	# }
	@pep_site;
}
