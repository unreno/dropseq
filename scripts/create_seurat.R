#!/usr/bin/env Rscript


#	if (!requireNamespace("BiocManager", quietly = TRUE))
#		install.packages("BiocManager")
#	BiocManager::install(c("devtools","Seurat","pryr","gdata","dplyr"),update = TRUE, ask = FALSE)


print("Loading libraries")

library(devtools)

library(Seurat)

library(pryr)
#	object_size(some_object)
#	humanReadable(mem_used())

library(gdata)


print(paste0("Starting at :",date(),":"))

print("Loading data from final.dge.txt.gz")

print(paste0("mem_used before read table :",humanReadable(mem_used()),": at :", date(),":"))

# load data
ds.data <- read.table("final.dge.txt.gz",row.names=1,header=T)
#	For 1, ( 17153128 Dec 30 01:52 final.dge.txt.gz )
#	Took over 1.5 hours and 76GB memory so far

print(paste0("mem_used after read table :",humanReadable(mem_used()),": at :", date(),":"))

print(paste0("nrow(ds.data) :",nrow(ds.data),":"))
print(paste0("ncol(ds.data) :",ncol(ds.data),":"))
print(paste0("object_size(ds.data) :",object_size(ds.data),":"))

print(paste0("mem_used before CreateSeurat :",humanReadable(mem_used()),": at :", date(),":"))
print("CreateSeuratObject")
#ds <- CreateSeuratObject(raw.data = ds.data, min.cells = 3,  min.genes = 200)	#	20180518

#	Substantial difference in Seurat version 2 and 3
#	counts instead of raw.data. same?
#	no min.genes. has min.features. same?
ds <- CreateSeuratObject( counts = ds.data, min.cells = 3, min.features = 200 )
print(ds)
print(paste0("mem_used after CreateSeurat, before rm(ds.data) :",humanReadable(mem_used()),": at :", date(),":"))


print("Removing raw ds.data")
rm(ds.data)

print(paste0("mem_used after rm(ds.data) :",humanReadable(mem_used()),": at :", date(),":"))

print(paste0("object_size(ds) :",object_size(ds),":"))

print("Saving ds")
save(ds, file="InitialSeuratObjectSample.RData")

print(paste0("Ending R script :",humanReadable(mem_used()),": at :", date(),":"))

