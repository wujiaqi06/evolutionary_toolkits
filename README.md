# evolutionary_toolkits
Tools to handle sequence data

## 1, Remove Gap Perl Suits
Remove gap regins from alignments that are placed in an input folder.
Alignments should have a name like "*.fas".


CutSeq_byReferenceGenome.pl   #Remove gaps that do not exist in the reference sequence.

RemoveGap_nuc.pl #Remove gaps from nucleotide sequence alignments.

RemoveGap_cds.pl #Remove gaps from codon alignments.

RemoveGap_pep.pl #Remove gaps from protein alignments.

Usage:
```
perl CutSeq_byReferenceGenome.pl -i Genome_folder -o output_folder -r reference_seq_name
perl RemoveGap_nuc.pl -i Genome_folder -o output_folder -d deleting_rate -t fasta/phy
perl RemoveGap_cds.pl -i Genome_folder -o output_folder -d deleting_rate -t fasta/phy
perl RemoveGap_pep.pl -i Genome_folder -o output_folder -d deleting_rate -t fasta/phy
```
