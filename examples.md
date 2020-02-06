# HEADER 1
## HEADER 2
### HEADER 3
#### HEADER 4
##### HEADER 5

+ List 
    - list
    - list 
    - list



## HORIZONTAL RULE
Three or more...

---

Hyphens (---)

*** 

Asterisks (***)

___

Underscores (___)



## SCRIPTS

<pre class="prettyprint"><code class="language-sh" style="background-color:333333">
#!/bin/bash

#SBATCH --job-name=star # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=60
#SBATCH --mem=32000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --reservation=workshop
#SBATCH --account=workshop
#SBATCH --array=1-16
#SBATCH --output=slurmout/star_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/star_%A_%a.err # File to which STDERR will be written

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`
REF="References/star.overlap100.gencode.v31"

outpath='02-STAR_alignment'
[[ -d ${outpath} ]] || mkdir ${outpath}
[[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}

echo "SAMPLE: ${sample}"

module load star/2.7.0e

call="STAR
     --runThreadN 8 \
     --genomeDir $REF \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesCommand zcat \
     --readFilesIn 01-HTS_Preproc/${sample}/${sample}_R1.fastq.gz 01-HTS_Preproc/${sample}/${sample}_R2.fastq.gz \
     --quantMode GeneCounts \
     --outFileNamePrefix ${outpath}/${sample}/${sample}_ \
     ${outpath}/${sample}/${sample}-STAR.stdout 2> ${outpath}/${sample}/${sample}-STAR.stderr"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime

</code></pre>

<pre class="prettyprint"><code class="language-py" style="background-color:333333">
#!/bin/python
print 'one'
print 'two'

if x == 1:
    print 'one'

cond1 = True
cond2 = False
if cond1 and cond2:
    # do something
    
</code></pre>

<pre class="prettyprint"><code class="language-R" style="background-color:333333">
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]

top.table <- data.frame(top.table,anno[match(top.table$Gene,anno$Gene.stable.ID.version),],logcpm[match(top.table$Gene,rownames(logcpm)),])

write.table(top.table, file = "A.C_v_B.C.txt", row.names = F, sep = "\t", quote = F)
    
</code></pre>

## Embedding PDFs
<object data="https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="http://yoursite.com/the.pdf">Download PDF</a>.</p>
    </embed>
</object>



## R CODE
- run the `alter_rmd.py` script to change some of the R formatting.
```bash
python alter_rmd.py -i Intro2R.md
```
- will produce a file with the same name `Intro2R_fixed.md`

### BEFORE
```r
# assign number 150 to variable a.
a <- 150
a
```

```
## [1] 150
```

### AFTER
```r
# assign number 150 to variable a.
a <- 150
a
```
<div class='r_output'> [1] 150
</div>

### TABLES FROM R
<table class="table table-striped table-hover table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:center;"> Description </th>
   <th style="text-align:center;"> R_function </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> Mean </td>
   <td style="text-align:center;"> mean() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Standard deviation </td>
   <td style="text-align:center;"> sd() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Variance </td>
   <td style="text-align:center;"> var() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Minimum </td>
   <td style="text-align:center;"> min() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Maximum </td>
   <td style="text-align:center;"> max() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Median </td>
   <td style="text-align:center;"> median() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Range of values: minimum and maximum </td>
   <td style="text-align:center;"> range() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sample quantiles </td>
   <td style="text-align:center;"> quantile() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Generic function </td>
   <td style="text-align:center;"> summary() </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Interquartile range </td>
   <td style="text-align:center;"> IQR() </td>
  </tr>
</tbody>
</table>



## OUTPUT (VERY WIDE) (actually this is a file)

<div class="output">Taxon_Name	MeanBootstrapValue	MeanLengthMerged	PercentageAsPairs	Total
d__Bacteria	0.997	421	0.0	126506
d__Bacteria;p__Acetothermia;c__Acetothermia_genera_incertae_sedis;o__Acetothermia_genera_incertae_sedis;f__Acetothermia_genera_incertae_sedis;g__Acetothermia_genera_incertae_sedis	0.56	421	0.0	1
d__Bacteria;p__Acidobacteria	0.605	424	0.0	1483
d__Bacteria;p__Acidobacteria;c__Acidobacteria_Gp1	0.983	403	0.0	1049
d__Bacteria;p__Acidobacteria;c__Acidobacteria_Gp10;o__Gp10;f__Gp10;g__Gp10	0.98	427	0.0	8312
d__Bacteria;p__Acidobacteria;c__Acidobacteria_Gp11;o__Gp11;f__Gp11;g__Gp11	0.792	406	0.0	321
d__Bacteria;p__Acidobacteria;c__Acidobacteria_Gp12;o__Gp12;f__Gp12;g__Gp12	0.999	403	0.0	34
d__Bacteria;p__Acidobacteria;c__Acidobacteria_Gp13;o__Gp13;f__Gp13;g__Gp13	0.998	423	0.0	13
d__Bacteria;p__Acidobacteria;c__Acidobacteria_Gp15;o__Gp15;f__Gp15;g__Gp15	0.961	414	0.0	1356
</div>