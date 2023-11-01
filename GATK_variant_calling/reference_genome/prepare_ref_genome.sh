# Index database sequences in the FASTA format
bwa index CD630.fasta

# Index reference sequence in the FASTA format
samtools faidx CD630.fasta

# Creates a sequence dictionary for a reference sequence
java -jar ../picard.jar CreateSequenceDictionary R=CD630.fasta O=CD630.dict

