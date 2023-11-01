curr_folder=$PWD
command="wait"
for path in $curr_folder/../wgs/*; do  # note that path is the absolute directory
    # if not a directory, skip
    [ -d "${path}" ] || continue
    base_folder=`basename $path`

    # enter the folder
    cd $path

    # get bam file with marked duplicates
    bamfile=(*.marked_duplicates.bam)
    base_fastq=${bamfile%.marked_duplicates.bam}

    # index bam file
    samtools index $bamfile

    # run haplotypecaller to call SNPs and indels simultaneously via local de-novo assembly of haplotypes
    $curr_folder/../gatk-4.1.9.0/gatk --java-options "-Xms4g -Xmx196g" HaplotypeCaller \
      -R $curr_folder/../reference_genome/CD630.fasta \
      -I $bamfile \
      -O $base_fastq.unfiltered.vcf.gz \
      -ploidy 1 \
      -ERC GVCF &

    command="$command $!"
done
eval $command
