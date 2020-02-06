#!/bin/bash

start=`date +%s`
echo $HOSTNAME

outpath='02-Salmon_alignment'
[[ -d ${outpath} ]] || mkdir ${outpath}

REF="References/salmon_gencode.v31.index"
GTF="References/gencode.v31.primary_assembly.annotation.gtf"

for sample in `cat samples.txt`
do
  [[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}
  echo "SAMPLE: ${sample}"

  call="salmon quant \
        --threads 8 \
        --index ${REF} \
        --libType A \
        --validateMappings \
        --geneMap ${GTF} \
        --output $outdir/$SAMPLE \
        -1 01-HTS_Preproc/${sample}/${sample}_R1.fastq.gz \
        -2 01-HTS_Preproc/${sample}/${sample}_R2.fastq.gz"

  echo $call
  eval $call
done

end=`date +%s`
runtime=$((end-start))
echo $runtime
