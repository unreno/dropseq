#!/usr/bin/env Rscript


print("Loading libraries")

#install.packages("devtools")
library(devtools)

#install_github("satijalab/seurat", ref = "3bd092a")	#	No. Doesn't include CreateSeuratObject and likely others.
#install.packages("httpuv")
#install.packages("Seurat")
library(Seurat)

#	install.packages("pryr")
library(pryr)
#	object_size(some_object)
#	humanReadable(mem_used())

library(gdata)


print(paste0("Starting at :",date(),":"))

print("Loading data from error_detected.dge.txt.gz")

print(paste0("mem_used before read table :",humanReadable(mem_used()),": at :", date(),":"))

# load data
ds.data <- read.table("error_detected.dge.txt.gz",row.names=1,header=T)
#	For 1, ( 17153128 Dec 30 01:52 error_detected.dge.txt.gz )
#	Took over 1.5 hours and 76GB memory so far

print(paste0("mem_used after read table :",humanReadable(mem_used()),": at :", date(),":"))

print(paste0("nrow(ds.data) :",nrow(ds.data),":"))
print(paste0("ncol(ds.data) :",ncol(ds.data),":"))
print(paste0("object_size(ds.data) :",object_size(ds.data),":"))


#ds <- CreateSeuratObject(raw.data = ds.data, min.cells = 3,  min.genes = 200, is.expr=1)
#
#	Error in if (!nreplace) return(x) : missing value where TRUE/FALSE needed
#	Calls: CreateSeuratObject ... suppressMessages -> withCallingHandlers -> [<- -> [<-.data.frame
#	In addition: Warning message:
#	In sum(i, na.rm = TRUE) : integer overflow - use sum(as.numeric(.))
#	Execution halted
#
#	Maybe, we need to remove "is.expr=1"?
#
#	ds <- CreateSeuratObject(raw.data = ds.data, min.cells = 3,  min.genes = 200)
#
#	Or filter more data at DigitalExpression stage in dropseq (i.e. MIN_NUM_GENES_PER_CELL=2000 instead of 100)?
#	Others kept only 0.1% of their data for Seurat.
#
#	Trying removing is.expr=1

print(paste0("mem_used before CreateSeurat :",humanReadable(mem_used()),": at :", date(),":"))

print("CreateSeuratObject")
#	ds <- CreateSeuratObject(raw.data = ds.data, min.cells = 3,  min.genes = 200)
#	ds <- CreateSeuratObject(raw.data = ds.data, min.cells = 3,  min.genes = 10)	#	20180314
ds <- CreateSeuratObject(raw.data = ds.data, min.cells = 3,  min.genes = 200)	#	20180518

print(paste0("mem_used after CreateSeurat, before rm(ds.data) :",humanReadable(mem_used()),": at :", date(),":"))

print("Removing raw ds.data")
rm(ds.data)

print(paste0("mem_used after rm(ds.data) :",humanReadable(mem_used()),": at :", date(),":"))

print(paste0("object_size(ds) :",object_size(ds),":"))

print("Saving ds")
save(ds, file="InitialSeuratObjectSample.RData")

print(paste0("Ending R script :",humanReadable(mem_used()),": at :", date(),":"))

