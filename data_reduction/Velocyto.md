
Login to tadpole and navigate to your directory on the share space.

```bash
cd /share/workshop/adv_scrna/$USER

srun -t 1-00:00:00 -c 4 -n 1 --mem 16000 --partition production --account adv_scrna_workshop --reservation adv_scrna_workshop  --pty /bin/bash
```

#  RNA Velocity measurement using Velocyto

RNA velocity is the time derivative of the gene expression state,
 [(La Manno et al., 2018)](https://www.nature.com/articles/s41586-018-0414-6) allows for the inference of the dynamic patterns in scRNA-seq data sets, by looking at the abundance of unspliced and spliced mRNA RNA in each cell, and modelling using a system of ordinary differential equations.

From a quantification point of view, RNA velocity analysis requires the generation of two count matrices, representing the spliced/processed (exonic) and unspliced/unprocessed RNA (intronic).

RNA velocity is a high-dimensional vector that predicts the future state of a cell on a timescale of hours.


### Example dataset

We'll use the example dataset from Mapping Comparison section, sample 654 from the paper

[_"Single-Cell RNA-seq Reveals Profound Alterations in Mechanosensitive Dorsal Root Ganglion Neurons With Vitamin E Deficiency"_](https://pubmed.ncbi.nlm.nih.gov/31733517/)


```bash
cd /share/workshop/adv_scrna/$USER/scrnaseq_processing
mkdir cd /share/workshop/adv_scrna/$USER/scrnaseq_processing/velocity
cd /share/workshop/adv_scrna/$USER/scrnaseq_processing/velocity
ln -s /share/biocore/workshops/2020_scRNAseq/velocity/654 .
```

we also need the gtf file used in creating the reference

```bash
ln -s /share/biocore/workshops/2020_scRNAseq/Reference/Mus_musculus.GRCm38.100.filtered.gtf .
```

### Next lets install the software [velocyto](https://velocyto.org/)

To install velocyto (a python application) we are going to use conda and a virtual environment

```bash
module load anaconda3
conda create -p velocyto
conda activate /share/workshop/adv_scrna/$USER/scrnaseq_processing/velocity/velocyto
```

If the environment 'activated' properly, than your prompt should look something like this.

```
(/share/workshop/adv_scrna/msettles/scrnaseq_processing/velocity/velocyto) msettles@tadpole:/share/workshop/adv_scrna/msettles/scrnaseq_processing/velocity$
```

If not you can try to initialize your bash
```bash
kinit
conda init bash
source /home/$USER/.bashrc
conda activate /share/workshop/adv_scrna/$USER/scrnaseq_processing/velocity/velocyto
```

Once your conda environment is properly activated you can then install the softare.
```bash
# install prerequisites
conda install numpy scipy cython numba matplotlib scikit-learn h5py click
# install velocyto
pip install velocyto
```

and then test the installation
```bash
velocyto --help
```

It is also recommended to use a repetative sequence mask file, this is easiest to obtain from the [UCSC genome browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Mouse&db=mm10&hgta_group=allTracks&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&position=chr12%3A56694976-56714605&hgta_outputType=primaryTable&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf
).

You want to download the mask in gtf format. For the sake of speed lets just link my copy. Unfortunately Ensembl doesn't provide a gtf file for th repeat sequences and we'd have to generate one ourselves which is beyond the scope of this workshop.

```bash
ln -s /share/biocore/workshops/2020_scRNAseq/velocity/mm10_rmsk.gtf .
```

###  Running Velocyto on Cellranger output

Lets first take a look at the help doc for run10x

```
velocyto run10x --help
```

Your dirctory should look something like This

```
drwxrwsr-x  3 msettles adv_scrna_workshop  6 May 28 07:18 .
drwxrwsr-x 13 msettles adv_scrna_workshop 21 May 28 06:39 ..
lrwxrwxrwx  1 msettles adv_scrna_workshop 51 May 28 07:07 654 -> /share/biocore/workshops/2020_scRNAseq/velocity/654
lrwxrwxrwx  1 msettles adv_scrna_workshop 76 May 28 07:10 Mus_musculus.GRCm38.100.gtf -> /share/biocore/workshops/2020_scRNAseq/Reference/Mus_musculus.GRCm38.100.filtered.gtf
drwxrwsr-x 22 msettles adv_scrna_workshop 22 May 28 06:58 velocyto
```

Ok, now we are ready to run this on our sample

```bash
cd /share/workshop/adv_scrna/$USER/scrnaseq_processing/velocity
velocyto run10x  654 Mus_musculus.GRCm38.100.filtered.gtf
```

### Output

Velocyto run10x simply produces a folder called velocyto in the sample directory with a single [loom](https://linnarssonlab.org/loompy/format/index.html) file in it, which contains the needed matrices for the analysis.

### More reading

A detailed study of the impact of quantification on RNA velocity estimates and interpretation, see [Soneson et al., 2020](https://www.biorxiv.org/content/10.1101/2020.03.13.990069v1).
