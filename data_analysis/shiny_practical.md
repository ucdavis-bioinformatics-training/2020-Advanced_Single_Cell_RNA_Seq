# App Practical Usage
This shiny app was created with the intention of working with biologists to extract meaning from the data by exploration.
Usage of the app is not for creating finalized "publish ready" images but rather a means for enabling an analysis of 
potential clustering based on the marker genes of interest (the biological question). The app enables easier interaction
with the data in order to explore what next steps can be taken towards producing "publish quality" images and which 
further analyses are relevant to the biological question. 

Manual intervention (such as clustering) with the Seurat objects is often needed following this exploration. The hope 
is for the app to enables some of these capabilities in the near future. Lets do some exercises to see how might use the
app to properly infer what custom clustering may need to be done.  



### Exercises
1) Use the single marker view to look at percent.mito at resolution 0.5 with a reduction of either tsne or umap. 
Which two groups cluster the most? What part of the app can you use to see if these two clusters are the most closely
related compared to the other clusters?

2) Extracting biological meaning in your data is an important part of scRNA-seq analysis. Typically this can be done by
looking through literature to confirm markers for cell types from your sample type. 
    -The data used for our analysis is from the follwoing paper: https://pubmed.ncbi.nlm.nih.gov/31733517/

    Check out the paper that this data is from, can you recognize any similar patterns? 
    
    Can you think of anything you might do to the Seurat object that might enable better analysis with the app based on the paper?

3) Run the app with the .RDS file that we created in the Anchoring Lecture. (Hint: you will need to change one line
in the app.R file). What are some differences you know when doing the comparison?



### A few extra things:
- Your dataset is huge? You have to run analysis on an HPC and visualizing/interacting with your data is tricky via CLI. 
AWS has huge instances you can run during the day for about 10$ with 100GB+ of RAM. You can just run it during working hours
when you can share the dataset with colleagues using a web hosted url, as it is billed hourly. 
    - If this is of interest to you feel free to reach out and check out the Shiny on AWS under Data Analysis in the docs.
    
- You are working with someone where it may likely be difficult for them to download packages, Rstudio, the app etc. and 
you need to collaborate without being able to both access the same computer. AWS can help with this as well!

- With both of the above statements in mind, we will be trying to allow hosting of the app ont he Genome Center HPC in the 
near future so keep an eye out for that!