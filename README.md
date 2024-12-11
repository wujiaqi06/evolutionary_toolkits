# evolutionary_toolkits
Tools to handle sequence data

##Remove Gap Perl Suits
CutSeq_byReferenceGenome.pl
RemoveGap_nuc.pl
RemoveGap_cds.pl
RemoveGap_pep.pl

Usage:
'''
perl CutSeq_byReferenceGenome.pl -i Genome_folder -o output_folder -r reference_seq_name
perl RemoveGap_nuc.pl -i Genome_folder -o output_folder -d deleting_rate -t fasta/phy
perl RemoveGap_pep.pl -i Genome_folder -o output_folder -d deleting_rate -t fasta/phy
perl RemoveGap_cds.pl -i Genome_folder -o output_folder -d deleting_rate -t fasta/phy
'''
