# Replace sample ID

cp ../wgs/cohort.filtered.gt.txt ../wgs/cohort.filtered.gt.final.txt
curr_folder=$PWD
for path in $curr_folder/../wgs/*; do  # note that path is the absolute directory
    # if not a directory, skip
    [ -d "${path}" ] || continue
    base_folder=`basename $path`

    # enter the folder
    cd $path

    # get sample names
    Reads1=(*R1.fastq.gz)
    Reads2=(*R2.fastq.gz)
    base_fastq=${Reads1%_R*}

    # get read group information
    header=$(zcat < $Reads1 | head -n 1)
    id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
    sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+\+[ATGCN]+$" | sed 's/+/_/g')

    # replace sample ID
    old_id=$id"_"$sm
    new_id=$base_fastq
    echo "$old_id, $new_id"
    sed -i '' "s/$old_id/$new_id/g" $curr_folder/../wgs/cohort.filtered.gt.final.txt

    cd $curr_folder
done
