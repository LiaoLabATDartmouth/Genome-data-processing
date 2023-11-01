# ./RepeatMasker/RepeatMasker -pa 2 -species bacteria -dir output ../reference_genome/CD630.fasta
# rmsk2bed < output/CD630.fasta.out > output/CD630.fasta.bed
python3 ./RepeatMasker/util/RM2Bed.py output/CD630.fasta.out output/CD630.fasta.bed
