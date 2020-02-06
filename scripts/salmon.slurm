#!/bin/bash

#SBATCH --array=1-16
#SBATCH --job-name=salmon # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=1440
#SBATCH --mem=20000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --reservation=workshop
#SBATCH --account=workshop
#SBATCH --output=slurmout/salmon_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/salmon_%A_%a.err # File to which STDERR will be written


start=`date +%s`
hostname

outdir="02-Salmon_alignment"
sampfile="samples.txt"
REF="References/salmon_gencode.v31.index"
GTF="References/gencode.v31.primary_assembly.annotation.gtf"

SAMPLE=`head -n ${SLURM_ARRAY_TASK_ID} $sampfile | tail -1`
R1="01-HTS_Preproc/$SAMPLE/${SAMPLE}_R1.fastq.gz"
R2="01-HTS_Preproc/$SAMPLE/${SAMPLE}_R2.fastq.gz"

echo $SAMPLE

if [ ! -e $outdir ]; then
    mkdir $outdir
fi

module load salmon

call="salmon quant \
      --threads 8 \
      --index ${REF} \
      --libType A \
      --validateMappings \
      --geneMap ${GTF} \
      --output $outdir/$SAMPLE \
      -1 $R1 \
      -2 $R2"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds
