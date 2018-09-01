#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied.", call.=FALSE)
}

last=args[length(args)]
rm(args)

#pattern=sprintf("^%s_",last)
#print(paste0("Will be searching pattern ", pattern))

print("Loading libraries")

library(devtools)
library(Seurat)
library(pryr)

date()

print("mem_used()")
mem_used()





#for( i in list.files(path="..", pattern=pattern) ){
#	print(paste0("Loading previously saved data from ",i,"."))
#	file=sprintf("../%s/InitialSeuratObjectSample.RData",i)
#	load( file )	#	loads "ds"
#
#	if( exists('merged') ){
#		merged <- MergeSeurat(object1 = merged, object2 = ds, add.cell.id2 = pattern)	#, project = "L001,L002")
#
#	} else {
#		merged = ds
#	}
#
#}

print(paste0("\n\nGrouping ",last,"\n\n"))


load(sprintf("../%s_L001/InitialSeuratObjectSample.RData",last))
print("Loaded L001")
ds1=ds
mem_used()
load(sprintf("../%s_L002/InitialSeuratObjectSample.RData",last))
print("Loaded L002")
ds2=ds
mem_used()
load(sprintf("../%s_L003/InitialSeuratObjectSample.RData",last))
print("Loaded L003")
ds3=ds
mem_used()
load(sprintf("../%s_L004/InitialSeuratObjectSample.RData",last))
print("Loaded L004")
ds4=ds
mem_used()
print("Removing ds")
rm(ds)
mem_used()
date()
print("Removing last")
rm(last)
mem_used()

ls(all.names = TRUE)

print("Merging L001 and L002")
ds1n2 <- MergeSeurat(object1 = ds1, object2 = ds2, add.cell.id1 = "L001", add.cell.id2 = "L002", project = "L001,L002", do.normalize = FALSE)
mem_used()
print("Removing L001")
rm(ds1)
mem_used()
print("Removing L002")
rm(ds2)
mem_used()
date()

print("Merging L003 and L004")
ds3n4 <- MergeSeurat(object1 = ds3, object2 = ds4, add.cell.id1 = "L003", add.cell.id2 = "L004", project = "L003,L004", do.normalize = FALSE)
mem_used()
print("Removing L003")
rm(ds3)
mem_used()
print("Removing L004")
rm(ds4)
mem_used()
date()

print("Merging L001/L002 and L003/L004")
ds <- MergeSeurat(object1 = ds1n2, object2 = ds3n4, project = "L001,L002,L003,L004", do.normalize = FALSE)
mem_used()
print("Removing merged L001/L002")
rm(ds1n2)
mem_used()
print("Removing merged L003/L004")
rm(ds3n4)
mem_used()
date()


ls(all.names = TRUE)


#stop("testing")

#	The following is the same as seurat.R without the initial read.table














print("Saving ds")
save(ds, file="InitialSeuratObjectSample.RData")
mem_used()
date()

print("Creating VlnPlot")
VlnPlot(object = ds, features.plot = c("nGene", "nUMI"), nCol = 2)
mem_used()
date()

print("Creating GenePlot")
GenePlot(object = ds, gene1 = "nUMI", gene2 = "nGene")
mem_used()
date()

print("NormalizeData")
ds <- NormalizeData(object = ds)
mem_used()
date()

print("ScaleDate")
ds <- ScaleData(object = ds)
mem_used()
date()

print("FindVariableGenes")
ds <- FindVariableGenes(object = ds, do.plot = TRUE)
mem_used()
date()

#	from http://satijalab.org/seurat/pbmc3k_tutorial.html

print("RunPCA")
ds <- RunPCA(object = ds, pc.genes = ds@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
mem_used()
date()

print("PrintPCA")
# Examine and visualize PCA results a few different ways
PrintPCA(object = ds, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
mem_used()
date()

print("VizPCA")
VizPCA(object = ds, pcs.use = 1:2)
mem_used()
date()

print("PCAPlot")
PCAPlot(object = ds, dim.1 = 1, dim.2 = 2)
mem_used()
date()

# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers
# that are strongly correlated with cellular heterogeneity, but may not have
# passed through variable gene selection.  The results of the projected PCA
# can be explored by setting use.full=T in the functions above
print("ProjectPCA")
ds <- ProjectPCA(object = ds, do.print = FALSE)
mem_used()
date()

print("Making a couple Heat Maps")

print("PCHeatmap")
PCHeatmap(object = ds, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
mem_used()
date()

print("PCHeatmap")
PCHeatmap(object = ds, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
mem_used()
date()

# NOTE: This process can take a long time for big datasets, comment out for
# expediency.  More approximate techniques such as those implemented in
# PCElbowPlot() can be used to reduce computation time
#ds <- JackStraw(object = ds, num.replicate = 100, do.print = FALSE)
#JackStrawPlot(object = ds, PCs = 1:12)
#PCElbowPlot(object = ds)





#	Warning messages:
#	1: In heatmap.2(data.use, Rowv = NA, Colv = NA, trace = "none", col = col.use,  :
#	  Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row dendogram.
#	2: In heatmap.2(data.use, Rowv = NA, Colv = NA, trace = "none", col = col.use,  :
#	  Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column dendogram.
#	3: In plot.window(...) : "dimTitle" is not a graphical parameter
#	4: In plot.xy(xy, type, ...) : "dimTitle" is not a graphical parameter
#	5: In title(...) : "dimTitle" is not a graphical parameter
#	Exception in thread "main" java.lang.NumberFormatException: For input string: "1e+05"
#		at java.lang.NumberFormatException.forInputString(NumberFormatException.java:65)
#		at java.lang.Integer.parseInt(Integer.java:580)
#		at java.lang.Integer.parseInt(Integer.java:615)
#		at ModularityOptimizer.readInputFile(ModularityOptimizer.java:184)
#		at ModularityOptimizer.main(ModularityOptimizer.java:74)
#	Error in file(file, "rt") : cannot open the connection
#	Calls: FindClusters -> RunModularityClustering -> read.table -> file
#	In addition: Warning message:
#	In file(file, "rt") :
#	  cannot open file '/Users/jakewendt/.R/output_50527.txt': No such file or directory
#	Execution halted







date()
print("Finding Clusters")

# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
print("FindClusters")
ds <- FindClusters(object = ds, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
mem_used()
date()

print("PrintFindClustersParams")
PrintFindClustersParams(object = ds)
mem_used()
date()

# While we do provide function-specific printing functions, the more general
# function to print calculation parameters is PrintCalcParams().


print("RunTSNE")
ds <- RunTSNE(object = ds, dims.use = 1:10, do.fast = TRUE, check_duplicates = FALSE)
mem_used()
date()
#Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  : 
#  Remove duplicates before running TSNE.

# note that you can set do.label=T to help label individual clusters

print("TSNEPlot")
TSNEPlot(object = ds)
mem_used()
date()



#	No point in this
#save(ds, file = "~/syryu/meddling.R")




#	Error in FindMarkers(object = ds, ident.1 = 1, min.pct = 0.25) : 
#	  No genes pass logfc.threshold threshold
#	Execution halted


date()
print("Finding Markers")

for (i in unique(FetchData(ds,"ident"))$ident){

	print(paste("Finding Markers for ident :",i,":", sep=""))
	# find all markers of cluster i
	date()
	cluster.markers <- FindMarkers(object = ds, ident.1 = i, min.pct = 0.25)
	print(x = head(x = cluster.markers, n = 10))

	write.table(x = head(x = cluster.markers, n = 10),
		file=paste("genelist_cluster.",i,".csv",sep=""),
		col.names=NA,
		sep="\t")

#		col.names=NA,	so first column has a header, albeit empty, as a placeholder. (if row.names = TRUE)

#	write.csv(x = head(x = cluster.markers, n = 10),
#		file=paste("genelist_cluster.",i,".csv",sep="") )
#
#		sep="\t",       NOT SETTABLE in write.csv. Need to use write.table
#		row.names=TRUE, DEFAULT
#		col.names=TRUE) UNKNOWN to write.csv (DEFAULT in write.table)

}



# find all markers distinguishing cluster 5 from clusters 0 and 3
#date()
#print("cluster5.markers")
#cluster5.markers <- FindMarkers(object = ds, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#print(x = head(x = cluster5.markers, n = 5))

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
#ds.markers <- FindAllMarkers(object = ds, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
#	Error in rbind(deparse.level, ...) : 
#	  numbers of columns of arguments do not match
#	In addition: Warning messages:
#	1: In data.frame(..., check.names = FALSE) :
#	  row names were found from a short variable and have been discarded
#	2: In data.frame(..., check.names = FALSE) :
#	  row names were found from a short variable and have been discarded
#	3: In data.frame(..., check.names = FALSE) :
#	  row names were found from a short variable and have been discarded


#	This can occur if individual clusters don't have any transcriptomic markers.
#	We'll improve in a future release, but you can either use FindMarkers, or lower the threshold for DE to solve.
#	Trying lower thresh.use
#ds.markers <- FindAllMarkers(object = ds, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.15)
#ds.markers <- FindAllMarkers(object = ds, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.05)
#ds.markers <- FindAllMarkers(object = ds, only.pos = TRUE, min.pct = 0.25 )
#ds.markers <- FindAllMarkers(object = ds, only.pos = TRUE, min.pct = 0.25, thresh.use = 0)
#ds.markers <- FindAllMarkers(object = ds, only.pos = TRUE, min.pct = 0.25, thresh.use = 10)
#
#	Nothing works on this data set


#print("Removing raw ds.data")
#rm(ds.data)

print("Saving all environment")
save(list=ls(all=TRUE), file="Sample.RData")



date()
print("mem_used()")
mem_used()
print("Ending R script")




#ds.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

#cluster1.markers <- FindMarkers(object = ds, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)






#	kinda specific to the example data so can't go any further until understand.

# you can plot raw UMI counts as well
#VlnPlot(object = ds, features.plot = c("MS4A1", "CD79A"))
#VlnPlot(object = ds, features.plot = c("NKG7", "PF4"), use.raw = TRUE, y.log = TRUE)
#FeaturePlot(object = ds, features.plot = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols.use = c("grey", "blue"), reduction.use = "tsne")
#
#top10 <- ds.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
#DoHeatmap(object = ds, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

