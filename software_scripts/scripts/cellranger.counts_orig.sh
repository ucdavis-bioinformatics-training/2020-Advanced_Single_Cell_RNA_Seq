#!/bin/bash

## Record the start time
start=`date +%s`

## Record the host being run on
echo "Hostname: $(eval hostname)"

## Where cellranger executable is located
## a) by loading a module
module load cellranger/3.1.0

## b) or, by placing the location of the executables on the path (edit to your location)
# export PATH=/share/pathtosoftware/bin/:$PATH

## c) or if they are already on the path, do nothing

## Set the parameters for the run
basepath=/share/workshop/adv_scrna/$USER/scrnaseq_processing
transcriptome="/share/workshop/adv_scrna/$USER/scrnaseq_processing/Reference/GRCm38.cellranger"
fastqs="${basepath}/00-RawData"

echo $basepath
echo $transcriptome
echo $fastq

for sample in `cat samples.txt`
do
    ## https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome
    ## Create the call
    call="cellranger count \
      --id=${sample} \
      --sample=${sample} \
      --transcriptome=${transcriptome} \
      --fastqs=${fastqs} \
      --localcores=4 \
      --localmem=15"

    ## Some other parameters that may be usefull/needed
    ## --expect-cells=NUM, number of cells expected
    ## --nosecondary, skip the unnecessary secondary analysis
    ## --r2-length=NUM, if your R2 qualities are really poor
    ## --chemistry=CHEM, should it fail chemistry detection

    ## Echo the call
    echo $call
    ## Evaluate the call
    eval $call
done

## Record the start time, and output runtime
end=`date +%s`
runtime=$((end-start))
echo $runtime
