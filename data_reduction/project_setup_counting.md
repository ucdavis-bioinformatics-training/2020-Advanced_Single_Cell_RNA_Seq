
The dataset used in this section is from a recently published mouse sample performed within the UCD DNA Technology and Bioinformtics Cores.  Its is a very small subset of reads from the experiment.

iScience. 2019 Nov 22;21:720-735. doi: 10.1016/j.isci.2019.10.064. Epub 2019 Oct 31.  
[_"Single-Cell RNA-seq Reveals Profound Alterations in Mechanosensitive Dorsal Root Ganglion Neurons With Vitamin E Deficiency"_](https://pubmed.ncbi.nlm.nih.gov/31733517/)  
Carrie J Finno, Janel Peterson, Mincheol Kang, Seojin Park, Matthew H Bordbari, Blythe Durbin-Johnson, Matthew Settles, Maria C Perez-Flores, Jeong H Lee, Ebenezer N Yamoah

# Single Cell Expression Data Setup

Let's set up a project directory for the analysis.

**1\.** First, create a directory for your user and the example project in the workshop directory:

```bash
cd
mkdir -p /share/workshop/adv_scrna/$USER/scrnaseq_processing
```

---

**2a\.** Next, go into that directory, create a raw data directory (we are going to call this 00-RawData) and cd into that directory. Lets then create symbolic links to the fastq files that contains the raw read data.

```bash
cd /share/workshop/adv_scrna/$USER/scrnaseq_processing
mkdir 00-RawData
cd 00-RawData/
ln -s /share/biocore/workshops/2020_scRNAseq/scExpression_data_small/* .
```

This directory now contains a folder for each "sample" (in this case just 1, 654_small) and the fastq files for each "sample" are in the sample folders.

**2b\.** lets create a sample sheet for the project, store sample names in a file called samples.txt

```bash
ls > ../samples.txt
cat ../samples.txt
```

---
**3a\.** Now, take a look at the raw data directory.

```bash
ls /share/workshop/adv_scrna/$USER/scrnaseq_processing/00-RawData
```

**3b\.** To see a list of the contents of each directory.

```bash
ls *
```

**3c\.** Lets get a better look at all the files in all of the directories.

```bash
ls -lah */*
```

---
**4a\.** View the contents of the files using the 'less' command, when gzipped used 'zless' (which is just the 'less' command for gzipped files, q to exit):

```bash
cd 654_small/
zless 654_small_S1_L008_I1_001.fastq.gz
zless 654_small_S1_L008_R1_001.fastq.gz
zless 654_small_S1_L008_R2_001.fastq.gz
```

Make sure you can identify which lines correspond to a single read and which lines are the header, sequence, and quality values. Press 'q' to exit less.

#### Questions
* What are the lengths of each read?
* Can you identify the 10X parts of the read, CellBC+UMI polyT, transcript?

**4b\.**

Then, let's figure out the number of reads in this file. A simple way to do that is to count the number of lines and divide by 4 (because the record of each read uses 4 lines). In order to do this use cat to output the uncompressed file and pipe that to "wc" to count the number of lines:


```bash
zcat 654_small_S1_L008_R1_001.fastq.gz | wc -l
```

Divide this number by 4 and you have the number of reads in this file.

#### Questions
* How many reads?
