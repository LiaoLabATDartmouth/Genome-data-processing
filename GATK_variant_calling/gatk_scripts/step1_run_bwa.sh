curr_folder=$PWD
for path in $curr_folder/../wgs/*; do  # note that path is the absolute directory
    # if not a directory, skip
    [ -d "${path}" ] || continue
    base_folder=`basename $path`

    # enter the folder
    cd $path

    # get pair-end read samples
    Reads1=(*R1.fastq.gz)
    Reads2=(*R2.fastq.gz)
    base_fastq=${Reads1%_R*} # remove the shortest path in Reads1 after _R*

    # get read group information
    header=$(zcat < $Reads1 | head -n 1)
    id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed "s/@//" | sed "s/:/_/g")
    sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+\+[ATGCN]+$" | sed "s/+/_/g")
    # echo "Read Group @RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA"

    # run bwa-mem to map short reads to a reference genome
    bwa mem \
    -M \
    -t 36 \
    -R "@RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA" \
    $curr_folder/../reference_genome/CD630.fasta \
    $Reads1 $Reads2 | gzip > $base_fastq.sam.gz

    # sort the outpu SAM file and output a BAM file
    # BAM files are sorted by reference coordinates
    samtools sort -O bam -T $base_fastq.sort -o $base_fastq.sort.bam -@ 36 $base_fastq.sam.gz

    # mark duplicates
    # locate and tag duplicate reads in a BAM file
    java -jar $curr_folder/../picard.jar MarkDuplicates I=$base_fastq.sort.bam O=$base_fastq.marked_duplicates.bam M=$base_fastq.marked_dup_metrics.txt

    # these intermediate files are not needed
    rm $base_fastq.sam.gz
    rm $base_fastq.sort.bam

    cd $curr_folder
done
