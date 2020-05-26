#!/bin/bash

## Record the start time
start=`date +%s`

## Record the host being run on
echo "Hostname: $(eval hostname)"

#export PATH=/share/biocore/workshops/2020_scRNAseq/salmon/bin/:$PATH

## Set the parameters for the run
basepath=/share/workshop/adv_scrna/${USER}/scrnaseq_processing
transcriptome=${basepath}/Reference/GRCm38.salmon_decoys
fastqs=${basepath}'/01-HTStream'
output=${basepath}'/654_small_hstream_salmon'
resources=${basepath}'/resources'

echo $basepath
echo $transcriptome
echo $fastq
echo $output
echo $resources

sample='654_small_htstream'
call="salmon alevin \
    -l ISR \
    -1 ${fastqs}/${sample}/${sample}_S1_L001_R1_001.fastq.gz \
    -2 ${fastqs}/${sample}/${sample}_S1_L001_R2_001.fastq.gz \
    --chromium  \
    -i ${transcriptome} \
    -p 4 \
    -o ${output} \
    --tgMap ${resources}/txp2gene.tsv"

echo ${call}
eval ${call}

## Record the start time, and output runtime
end=`date +%s`
runtime=$((end-start))
echo $runtime
