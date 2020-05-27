#!/bin/env bash

## Record the start time
start=`date +%s`

## Record the host being run on
echo "Hostname: $(eval hostname)"

# If not already on the path
#export PATH=/share/workshop/adv_scrna/$USER/HTStream/bin:$PATH

## Set the parameters for the run
basepath=/share/workshop/adv_scrna/$USER/scrnaseq_processing
resources=${basepath}'/resources'
fastqpath=${basepath}'/00-RawData/654_small'
outpath=${basepath}'/01-HTStream'

echo $basepath
echo $resources
echo $fastqpath
echo $output

for sample in `cat samples.txt`
do
    echo "SAMPLE: ${sample}"
    mkdir -p ${outpath}/${sample}

    call="hts_Stats -L ${outpath}/${sample}/${sample}_scRNA.log -N 'compute stats on original dataset' \
        -1 ${fastqpath}/${sample}_S*_R1_001.fastq.gz \
        -2 ${fastqpath}/${sample}_S*_R2_001.fastq.gz  | \
    hts_SeqScreener -A ${outpath}/${sample}/${sample}_scRNA.log -N 'screen for PhiX because I always do' \
        --check-read-2 | \
    hts_Overlapper -A ${outpath}/${sample}/${sample}_scRNA.log -N 'overlap reads' | \
    extract_BC-UMI.py --extract --read 1 --length  28 | \
    hts_PolyATTrim -A ${outpath}/${sample}/${sample}_scRNA.log -N 'trim 3 prime plolyA' \
        --skip_polyA  \
        --no-right -x 100 | \
    hts_NTrimmer -A ${outpath}/${sample}/${sample}_scRNA.log -N 'Remove any N characters' | \
    hts_QWindowTrim -A ${outpath}/${sample}/${sample}_scRNA.log -N 'Quality window trim' | \
    hts_LengthFilter -A ${outpath}/${sample}/${sample}_scRNA.log -N 'Removed any read shorter than 50bp' \
        -m 50 -s | \
    hts_SeqScreener -A ${outpath}/${sample}/${sample}_scRNA.log -N 'Screen out any potential adapter dimers' \
        --check-read-2 \
        -s ${resources}/screen.fa | \
    hts_SeqScreener -A ${outpath}/${sample}/${sample}_scRNA.log -N 'count the number of rRNA reads'\
        -r \
        --check-read-2 \
        -s ${resources}/mouse_rrna.fasta | \
    extract_BC-UMI.py --insert --read 1 | \
    hts_Stats -A ${outpath}/${sample}/${sample}_scRNA.log -N 'final stats' \
        -f ${outpath}/${sample}/${sample} -F"

    echo $call
    eval $call
done

end=`date +%s`
runtime=$((end-start))
echo $runtime
