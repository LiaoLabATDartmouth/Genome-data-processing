# Combine unfiltered variants from all samples
command="../gatk-4.1.9.0/gatk CombineGVCFs -R ../reference_genome/CD630.fasta"
curr_folder=$PWD
for path in $curr_folder/../wgs/*; do  # note that path is the absolute directory
    # if not a directory, skip
    [ -d "${path}" ] || continue
    base_folder=`basename $path`

    # enter the folder
    cd $path

    # get bam file with marked duplicates
    bamfile=(*.marked_duplicates.bam)
    base_fastq=${bamfile%.marked_duplicates.bam}
    command="$command --variant $curr_folder/../wgs/$base_folder/$base_fastq.unfiltered.vcf.gz"

    cd $curr_folder
done
command="$command -O ../wgs/cohort.unfiltered.vcf.gz"
eval $command

# Perform joint genotyping on all samples pre-called with HaplotypeCaller
$curr_folder/../gatk-4.1.9.0/gatk --java-options "-Xms4g -Xmx196g" GenotypeGVCFs \
    -R ../reference_genome/CD630.fasta \
    -V ../wgs/cohort.unfiltered.vcf.gz \
    -ploidy 1 \
    -O ../wgs/cohort.unfiltered.gt.vcf.gz
gunzip -k ../wgs/cohort.unfiltered.gt.vcf.gz
