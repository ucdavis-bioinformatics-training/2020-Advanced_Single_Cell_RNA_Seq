#!/bin/bash

## assumes star version 2.7.0e

start=`date +%s`
echo $HOSTNAME

outpath="References"
[[ -d ${outpath} ]] || mkdir ${outpath}

cd ${outpath}
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.pc_transcripts.fa.gz
gunzip gencode.v31.pc_transcripts.fa.gz
PC_FASTA="gencode.v31.pc_transcripts.fa"
INDEX="salmon_gencode.v31.index"

module load salmon
call="salmon index -i ${INDEX} -k 31 --gencode -p 8 -t ${PC_FASTA}"
echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
