#!/bin/bash

#SBATCH --job-name=htstream # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --time=60
#SBATCH --mem=3000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --reservation=workshop
#SBATCH --account=workshop
#SBATCH --array=1-16
#SBATCH --output=slurmout/htstream_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/htstream_%A_%a.err # File to which STDERR will be written

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`

outpath='01-HTS_Preproc'
[[ -d ${outpath} ]] || mkdir ${outpath}
[[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}

echo "SAMPLE: ${sample}"

module load htstream/1.1.0

call="hts_Stats -L ${outpath}/${sample}/${sample}_htsStats.log -1 00-RawData/${sample}/*R1* -2 00-RawData/${sample}/*R2* | \
      hts_SeqScreener -A -L ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_SeqScreener -s References/human_rrna.fasta -r -A -L ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_SuperDeduper -e 250000 -A -L ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_AdapterTrimmer -n -A -L ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_QWindowTrim -n -A -L ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_NTrimmer -n -A -L ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_CutTrim -n -m 50 -A -L ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_Stats -A -L ${outpath}/${sample}/${sample}_htsStats.log -f ${outpath}/${sample}/${sample}"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
