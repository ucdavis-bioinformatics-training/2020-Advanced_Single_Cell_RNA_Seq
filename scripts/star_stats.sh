#!/bin/bash

echo -en "sample_names" > names.txt
echo -en "total_in_feature\t" > totals.txt
cat 02-STAR_alignment/*/*ReadsPerGene.out.tab | head -4 | cut -f1 > stats.txt
cat samples.txt | while read sample; do
    echo ${sample}
    echo -en "\t${sample}" >> names.txt
    head -4 02-STAR_alignment/${sample}/${sample}_ReadsPerGene.out.tab | cut -f4 > temp1
    paste stats.txt temp1 > temp2
    mv temp2 stats.txt
    tail -n +5 02-STAR_alignment/${sample}/${sample}_ReadsPerGene.out.tab | cut -f4 | \
        perl -ne '$tot+=$_ }{ print "$tot\t"' >> totals.txt
done
echo -en "\n" >> names.txt
cat names.txt stats.txt totals.txt > temp1
mv temp1 summary_star_alignments.txt
rm stats.txt
rm names.txt
rm totals.txt
