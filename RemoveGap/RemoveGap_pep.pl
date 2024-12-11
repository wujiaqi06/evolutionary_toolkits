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
getopts('i:o:d:t:', \%opts);
my $input_folder = $opts{'i'} or die "use: perl $0 -i Genome_folder -o output_folder -d deleting_rate -t fasta\/phy\n";
my $output_folder = $opts{'o'} or die "use: perl $0 -i Genome_folder -o output_folder -d deleting_rate -t fasta\/phy\n";
my $del_missing_value = $opts{'d'} or die "use: perl $0 -i Genome_folder -o output_folder -d deleting_rate -t fasta\/phy\n";
my $output_type = $opts{'t'} or die "use: perl $0 -i Genome_folder -o output_folder -d deleting_rate -t fasta\/phy\n";

#my $output_folder = "DelAmbSites";
#my $input_folder = "test";

mkdir "$output_folder" or die "Can not create folder $output_folder:$!\n";
open DELIND, ">./$output_folder/Delete_info.txt";
chdir "./$input_folder";
my @all_seq_file = glob '*.fas*';

foreach my $seq (@all_seq_file){
	my %read_sequences = &read_fasta("$seq");

	### generation the missing value tables: "-" as 0, no-missing as "1", saving in @missing
	my @sequences = values %read_sequences;
	my @sequence_name = sort keys %read_sequences;
	my $sequence_number = $#sequence_name +1;
	my $sequence_length = length($sequences[0]);
	my $missing_value = 0 x $sequence_length;
	my @missing = split(//, $missing_value);
	#print "missing is @missing\n";
	foreach my $i (@sequences) {
		my $temp_missing = "";
		my @seq_site = split(//,$i);
		#say @seq_site;
		foreach my $j (@seq_site){
			if (($j eq "-") | ($j eq "X")){
				$temp_missing .= 0;
			} else {$temp_missing .= 1;
				}
		}
		my @temp_array = split(//, $temp_missing);
		foreach my $j (0..$#missing){
			$missing[$j] = $missing[$j] + $temp_array[$j];
		}
		#print "temp_missing is @temp_array\n";
	}

	#print "missing_value is @missing\n";

	##getting the index of sites which should be deleted. saving in @del_index
	#my $del_missing_value = 0.8;
	my $cut_number = int (($#sequence_name+1) * $del_missing_value);
	#print "$cut_number\n";
	my @del_index;
	foreach my $i (0..$#missing){
		if ($missing[$i] < $cut_number){
			push @del_index, $i;
		}
	}
	my $delete_length = $#del_index + 1;
	my $remaining = $sequence_length - $delete_length;
	print "Delete $delete_length sites, remaining $remaining sites\n\n";
	print DELIND "$seq: original length: $sequence_length bp, deleted $delete_length bp, remaining $remaining bp.\n";

	## Delete the sites written in @del_index file
	@del_index = reverse (@del_index);
	if ($output_type =~ /fas(ta)?/){
		open OUTPUT, ">../$output_folder/$seq.fas\n";
		foreach my $seq_name (@sequence_name){
			my @site = split(//, $read_sequences{$seq_name});
			#my @site_del = splice @site, @del_index;
			foreach my $del (@del_index){
				splice(@site, $del, 1);
			}
			#print "length of site is $#site\n";
			#my $new_seq = join(/:/,@site);
			print OUTPUT ">$seq_name\n";
			print OUTPUT @site;
			print OUTPUT "\n";
		}
		close OUTPUT;
	} elsif($output_type eq "phy"){
		open OUTPUT, ">../$output_folder/$seq.phy\n";
		print OUTPUT "$sequence_number $remaining\n";
		foreach my $seq_name (@sequence_name){
			my @site = split(//, $read_sequences{$seq_name});
			#my @site_del = splice @site, @del_index;
			foreach my $del (@del_index){
				splice(@site, $del, 1);
			}
			#print "length of site is $#site\n";
			#my $new_seq = join(/:/,@site);
			print OUTPUT "$seq_name          ";
			print OUTPUT @site;
			print OUTPUT "\n";
		}
		close OUTPUT;
	}

}
close DELIND;
#print "del_index is @del_index\n";
#print "missing is $missing[1169]\n";



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


