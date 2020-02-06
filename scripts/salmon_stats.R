## first need to install the R packages JSON

#lets install to current working directory and then load
suppressWarnings(dir.create("r_lib"))
new_rlib = file.path(getwd(),"r_lib")

if (!any(rownames(installed.packages()) == "rjson")) install.packages("rjson", repos='http://cran.us.r-project.org', lib=new_rlib)

require(rjson, lib.loc=new_rlib)

samples <- readLines("samples.txt")
salmon_dir <- "02-Salmon_alignment"

salmon_json <- lapply(samples, function(s) {
  suppressWarnings(fromJSON(paste(readLines(file.path(salmon_dir,s,"aux_info/meta_info.json")),collapse="")))
})

if (length(salmon_json) != length(samples)) stop("not all samples have json log files")

## grap stats
num_processed <- sapply(salmon_json, function(x) x$"num_processed")
num_mapped <- sapply(salmon_json, function(x) x$"num_mapped")
num_decoy_fragments <- sapply(salmon_json, function(x) x$"num_decoy_fragments")
num_fragments_filtered_vm <- sapply(salmon_json, function(x) x$"num_fragments_filtered_vm")
num_alignments_below_threshold <- sapply(salmon_json, function(x) x$"num_alignments_below_threshold_for_mapped_fragments_vm")

LongTable <- data.frame(
    sample_name=samples,
    Num_processed=num_processed,
    Num_mapped=num_mapped,
    Num_decoy_fragments=num_decoy_fragments,
    Num_fragments_filtered=num_fragments_filtered_vm,
    Num_alignments_below_threshold=num_alignments_below_threshold,
    Percent_mapped=round((num_mapped)/(num_processed)*100,2)
)

write.table(t(LongTable),"summary_salmon_alignments.txt",sep="\t",row.names=T,col.names=F,quote=F)
