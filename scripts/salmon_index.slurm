#!/bin/bash

#SBATCH --job-name=salmon_index # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=60
#SBATCH --mem=15000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --reservation=workshop
#SBATCH --account=workshop
#SBATCH --output=slurmout/salmon-index_%A.out # File to which STDOUT will be written
#SBATCH --error=slurmout/salmon-index_%A.err # File to which STDERR will be written

start=`date +%s`
echo $HOSTNAME

outpath="References"

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
