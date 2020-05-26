#!/bin/bash

## Record the start time
start=`date +%s`

## Record the host being run on
echo "Hostname: $(eval hostname)"

## Where cellranger and bcl2fastq executables are located
## a) by loading a module
module load star/2.7.3a

## b) or, by placing the location of the executables on the path (edit to your location)
# export PATH=/share/pathtosoftware/bin/:$PATH

## c) or if they are already on the path, do nothing

## Set the parameters for the run
basepath='/share/workshop/adv_scrna/$USER/scrnaseq_processing'
transcriptome=${basepath}'/Reference/GRCm38.star'
fastqs=${basepath}'/01-HTStream'
output=${basepath}'/654_small_hstream_star'
resources=${basepath}'/resources'

echo $basepath
echo $transcriptome
echo $fastq
echo $output
echo $resources

sample='654_small_htstream'

call="STAR --runThreadN 4 \
    --runDirPerm All_RWX \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist ${resources}/737K-august-2016.txt \
    --soloCBstart 1 \
    --soloCBlen 16 \
    --soloUMIstart 17 \
    --soloUMIlen 12 \
    --soloBarcodeReadLength 0 \
    --soloCBmatchWLtype 1MM_multi \
    --soloUMIdedup 1MM_All \
    --soloStrand Forward \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 20 \
    --genomeDir ${transcriptome} \
    --outFileNamePrefix ${output}/star_654_preproc_ \
    --outSAMtype BAM SortedByCoordinate \
    --soloFeatures Gene Velocyto \
    --readFilesCommand zcat \
    --readFilesIn ${fastq}/${sample}/${sample}_S1_L008_R2_001.fastq.gz ${fastq}/${sample}/${sample}/654_preproc_S1_L008_R1_001.fastq.gz"

echo ${call}
eval ${call}

## Record the start time, and output runtime
end=`date +%s`
runtime=$((end-start))
echo $runtime
