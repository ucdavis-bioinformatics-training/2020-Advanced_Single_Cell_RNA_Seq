
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

We'll use the example dataset from Mapping Comparison section from the paper, sample 654_small

[_"Single-Cell RNA-seq Reveals Profound Alterations in Mechanosensitive Dorsal Root Ganglion Neurons With Vitamin E Deficiency"_](https://pubmed.ncbi.nlm.nih.gov/31733517/)


and the cellranger count results folder from the [Mapping](scMapping) section.

```
/share/workshop/adv_scrna/msettles/scrnaseq_processing/654_small
```

```bash
cd /share/workshop/adv_scrna/$USER/scrnaseq_processing
```

we also need the gtf file used in creating the reference, which should be here:

```
/share/workshop/adv_scrna/msettles/scrnaseq_processing/Reference/Mus_musculus.GRCm38.100.filtered.gtf .
```

### Next lets install the software [velocyto](https://velocyto.org/)

To install velocyto (a python application) we are going to use conda and a virtual environment

```bash
cd /share/workshop/adv_scrna/$USER
module load anaconda3
conda create -p velocyto
conda activate /share/workshop/adv_scrna/$USER/velocyto
```

If the environment 'activated' properly, than your prompt should look something like this.

```
(/share/workshop/adv_scrna/msettles/velocyto) msettles@tadpole:/share/workshop/adv_scrna/msettles$
```

If not you can try to initialize your bash
```bash
kinit
conda init bash
source /home/$USER/.bashrc
conda activate /share/workshop/adv_scrna/$USER/velocyto
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

You want to download the mask in gtf format. Unfortunately, Ensembl doesn't provide a gtf file for th repeat sequences and we'd have to generate one ourselves (write a script to do so) which is beyond the scope of this workshop.

###  Running Velocyto on Cellranger output

Lets first take a look at the help doc for run10x

```
velocyto run10x --help
```

Ok, now we are ready to run this on our sample

```bash
cd /share/workshop/adv_scrna/$USER/scrnaseq_processing
velocyto run10x  654_small References/Mus_musculus.GRCm38.100.filtered.gtf
```

### Output

Velocyto run10x simply produces a folder called velocyto in the sample directory with a single [loom](https://linnarssonlab.org/loompy/format/index.html) file in it, which contains the needed matrices for the analysis.

The output folder 654_small, now has a new folder called velocyto
```
drwxrwsr-x 5 msettles biocore   20 May 28 07:22 .
drwxrwsr-x 3 msettles biocore    4 May 28 07:17 ..
-rw-rw-r-- 1 msettles biocore 4.2M May 28 06:29 654.mri.tgz
-rw-rw-r-- 1 msettles biocore  287 May 28 06:29 _cmdline
-rw-rw-r-- 1 msettles biocore  66K May 28 06:29 _filelist
-rw-r--r-- 1 msettles biocore 969K May 28 06:29 _finalstate
-rw-r--r-- 1 msettles biocore  719 May 28 06:29 _invocation
-rw-r--r-- 1 msettles biocore    5 May 28 06:29 _jobmode
-rw-r--r-- 1 msettles biocore  69K May 28 06:27 _log
-rw-r--r-- 1 msettles biocore  66K May 28 06:29 _mrosource
drwxrwsr-x 5 msettles biocore   14 May 28 06:28 outs
-rw-r--r-- 1 msettles biocore 577K May 28 06:29 _perf
drwxrwsr-x 6 msettles biocore    6 May 28 06:26 SC_RNA_COUNTER_CS
-rw-rw-r-- 1 msettles biocore  17K May 28 06:29 _sitecheck
-rw-r--r-- 1 msettles biocore    2 May 28 06:29 _tags
-rw-r--r-- 1 msettles biocore   51 May 28 06:27 _timestamp
-rw-r--r-- 1 msettles biocore   36 May 28 06:29 _uuid
-rw-r--r-- 1 msettles biocore 132K May 28 06:27 _vdrkill
drwxrwsr-x 2 msettles biocore    3 May 28 08:50 velocyto    <-----
-rw-r--r-- 1 msettles biocore   61 May 28 06:27 _versions
```

and inside that folder is the loom file.

```
drwxrwsr-x 2 msettles biocore   3 May 28 08:50 .
drwxrwsr-x 5 msettles biocore  20 May 28 07:22 ..
-rw-rw-r-- 1 msettles biocore 56M May 28 08:52 654.loom
```

I expect that the 654_small dataset will fail as their aren't enough reads, but on a full dataset it produces the loom file.

### More reading

A detailed study of the impact of quantification on RNA velocity estimates and interpretation, see [Soneson et al., 2020](https://www.biorxiv.org/content/10.1101/2020.03.13.990069v1).
