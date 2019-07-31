#!/usr/bin/env Rscript

#	if (!requireNamespace("BiocManager", quietly = TRUE))
#		install.packages("BiocManager")
#	BiocManager::install(c("devtools","Seurat","pryr","gdata","dplyr","umap"),update = TRUE, ask = FALSE)

print("Loading libraries")

library(devtools)

library(Seurat)

library(pryr)
#	object_size(some_object)
#	humanReadable(mem_used())

library(gdata)

library(dplyr)
#	for function "%>%"


print(paste0("Starting at :",date(),":"))

print("Loading InitialSeuratObjectSample.RData")
print(paste0("mem_used before loading seurat object :",humanReadable(mem_used()),": at :", date(),":"))
load("InitialSeuratObjectSample.RData")
print(paste0("mem_used after loading seurat object :",humanReadable(mem_used()),": at :", date(),":"))

#	from http://satijalab.org/seurat/Seurat_AlignmentTutorial.html
#
#	then https://satijalab.org/seurat/essential_commands.html
#	 and https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html





print("PercentageFeatureSet")
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
ds[["percent.mt"]] <- PercentageFeatureSet(ds, pattern = "^MT-")


print("VlnPlot")
# Visualize QC metrics as a violin plot
VlnPlot(ds, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


print("CombinePlots")
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(ds, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ds, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


ds <- subset(ds, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


print("NormalizeData")
ds <- NormalizeData(ds, normalization.method = "LogNormalize", scale.factor = 10000)
#	same as just
#	ds <- NormalizeData(ds)


print("FindVariableFeatures")
ds <- FindVariableFeatures(ds, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ds), 10)


print("CombinePlots")
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ds)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))


print("ScaleData")
all.genes <- rownames(ds)
ds <- ScaleData(ds, features = all.genes)


print("RunPCA")
ds <- RunPCA(ds, features = VariableFeatures(object = ds))

# Examine and visualize PCA results a few different ways
print(ds[["pca"]], dims = 1:5, nfeatures = 5)


print("VizDimLoadings")
VizDimLoadings(ds, dims = 1:2, reduction = "pca")


print("DimPlot")
DimPlot(ds, reduction = "pca")


print("DimHeatmap")
DimHeatmap(ds, dims = 1, cells = 500, balanced = TRUE)


print("DimHeatmap")
DimHeatmap(ds, dims = 1:15, cells = 500, balanced = TRUE)


print("JackStraw")
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
ds <- JackStraw(ds, num.replicate = 100)
ds <- ScoreJackStraw(ds, dims = 1:20)


print("JackStrawPlot")
JackStrawPlot(ds, dims = 1:15)


print("ElbowPlot")
ElbowPlot(ds)


print("FindNeighbors")
ds <- FindNeighbors(ds, dims = 1:10)
ds <- FindClusters(ds, resolution = 0.5)


# Look at cluster IDs of the first 5 cells
head(Idents(ds), 5)


print("RunUMAP")
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
ds <- RunUMAP(ds, dims = 1:10)


print("DimPlot")
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(ds, reduction = "umap")


saveRDS(ds, file = "pbmc_tutorial.rds")


print("FindMarkers")
# find all markers of cluster 1
cluster1.markers <- FindMarkers(ds, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

#	Only works if 5 clusters
#	print("FindMarkers")
#	# find all markers distinguishing cluster 5 from clusters 0 and 3
#	cluster5.markers <- FindMarkers(ds, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#	head(cluster5.markers, n = 5)


print("FindAllMarkers")
# find markers for every cluster compared to all remaining cells, report only the positive ones
ds.markers <- FindAllMarkers(ds, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ds.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


print("FindMarkers")
cluster1.markers <- FindMarkers(ds, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#
#	This is very specific to the data.
#	  None of the requested variables were found: MS4A1, CD79A
#	VlnPlot(ds, features = c("MS4A1", "CD79A"))
#
#	This too.
# you can plot raw counts as well
#	  None of the requested variables were found: NKG7, PF4
#	VlnPlot(ds, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
#
#	Same here
#	Error: None of the requested features were found: MS4A1, GNLY, CD3E, CD14, FCER1A, FCGR3A, LYZ, PPBP, CD8A in slot data
#	FeaturePlot(ds, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
#

print("top10 DoHeatmap")
top10 <- ds.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(ds, features = top10$gene) + NoLegend()

#
#	Seems specific too
#new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
#names(new.cluster.ids) <- levels(ds)
#ds <- RenameIdents(ds, new.cluster.ids)
#DimPlot(ds, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#


saveRDS(ds, file = "pbmc3k_final.rds")

##################################################









print("NormalizeData")
ds <- NormalizeData(object = ds)
print(paste0("mem_used after NormalizeData :",humanReadable(mem_used()),": at :", date(),":"))



#	Not in Seurat version 3
#print("FindVariableGenes")
##ds <- FindVariableGenes(object = ds, do.plot = FALSE)
#ds <- FindVariableGenes(object = ds, do.plot = TRUE)
#print(paste0("mem_used after FindVariableGenes :",humanReadable(mem_used()),": at :", date(),":"))
#	But this is
print("FindVariableFeatures")	#	must be run before calling RunPCA
ds <- FindVariableFeatures(object = ds, do.plot = TRUE)
print(paste0("mem_used after FindVariableFeatures :",humanReadable(mem_used()),": at :", date(),":"))


print("ScaleData")
ds <- ScaleData(object = ds)
print(paste0("mem_used after ScaleData :",humanReadable(mem_used()),": at :", date(),":"))



#	print("RunPCA")
#	ds <- RunPCA(object = ds)
#	print("FindNeighbors")
#	ds <- FindNeighbors(object = ds)
#	print("FindClusters")
#	ds <- FindClusters(object = ds)
#	print("RunTSNE")
#	ds <- RunTSNE(object = ds)
#	print("DimPlot")
#	DimPlot(object = ds, reduction = "tsne")



#	https://satijalab.org/seurat/essential_commands.html
#ds.counts <- Read10X(data.dir = "~/Downloads/pbmc3k/filtered_gene_bc_matrices/hg19/")
#ds <- CreateSeuratObject(counts = ds.counts)
#ds <- NormalizeData(object = ds)
#ds <- FindVariableFeatures(object = ds)
#ds <- ScaleData(object = ds)
#ds <- RunPCA(object = ds)
#ds <- FindNeighbors(object = ds)
#ds <- FindClusters(object = ds)
#ds <- RunTSNE(object = ds)
#DimPlot(object = ds, reduction = "tsne")










#	from http://satijalab.org/seurat/pbmc3k_tutorial.html

print("RunPCA")
ds <- RunPCA(object = ds, pc.genes = ds@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
print(paste0("mem_used after RunPCA :",humanReadable(mem_used()),": at :", date(),":"))

print("PrintPCA")
# Examine and visualize PCA results a few different ways
print(ds[["pca"]], dims = 1:5, nfeatures = 5)
#
#	Not in Version 3
#
#PrintPCA(object = ds, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
#print(paste0("mem_used after PrintPCA :",humanReadable(mem_used()),": at :", date(),":"))

#print("VizPCA")
#
#	Not in Version 3
#
#VizPCA(object = ds, pcs.use = 1:2)
#print(paste0("mem_used after VizPCA :",humanReadable(mem_used()),": at :", date(),":"))
#	print("VizDimLoadings")
#	VizDimLoadings(ds, dims = 1:2, reduction = "pca")

print("DimPlot")
DimPlot(ds, reduction = "pca")
#	save as
#	print("PCAPlot")
#	PCAPlot(object = ds, dim.1 = 1, dim.2 = 2)
#	print(paste0("mem_used after PCAPlot :",humanReadable(mem_used()),": at :", date(),":"))

# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers
# that are strongly correlated with cellular heterogeneity, but may not have
# passed through variable gene selection.  The results of the projected PCA
# can be explored by setting use.full=T in the functions above
#
#	Not in Version 3
#
#print("ProjectPCA")
#ds <- ProjectPCA(object = ds, do.print = FALSE)
#print(paste0("mem_used after ProjectPCA :",humanReadable(mem_used()),": at :", date(),":"))



#	Use those above
#	print("Making a couple Heat Maps")
#	
#	#print("PCHeatmap")
#	#PCHeatmap(object = ds, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
#	#print(paste0("mem_used after PCHeatmap :",humanReadable(mem_used()),": at :", date(),":"))
#	print("DimHeatmap")
#	DimHeatmap(ds, dims = 1, cells = 500, balanced = TRUE)
#	
#	#print("PCHeatmap")
#	#PCHeatmap(object = ds, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
#	#print(paste0("mem_used after PCHeatmap :",humanReadable(mem_used()),": at :", date(),":"))
#	print("DimHeatmap")
#	DimHeatmap(ds, dims = 1:15, cells = 500, balanced = TRUE)



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
#ds <- FindClusters(object = ds, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
#print(paste0("mem_used after FindClusters :",humanReadable(mem_used()),": at :", date(),":"))
ds <- FindNeighbors(ds, dims = 1:10)
ds <- FindClusters(ds, resolution = 0.5)



#
#	Not in Version 3
#
#print("PrintFindClustersParams")
#PrintFindClustersParams(object = ds)
#print(paste0("mem_used after PrintFindClustersParams :",humanReadable(mem_used()),": at :", date(),":"))

# While we do provide function-specific printing functions, the more general
# function to print calculation parameters is PrintCalcParams().

print("RunTSNE")
ds <- RunTSNE(object = ds, dims.use = 1:10, do.fast = TRUE, check_duplicates = FALSE)
print(paste0("mem_used after RunTSNE :",humanReadable(mem_used()),": at :", date(),":"))

#Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  :
#  Remove duplicates before running TSNE.

# note that you can set do.label=T to help label individual clusters

print("TSNEPlot")
TSNEPlot(object = ds)
print(paste0("mem_used after TSNEPlot :",humanReadable(mem_used()),": at :", date(),":"))





date()
print("Finding Markers")

for (i in unique(FetchData(ds,"ident"))$ident){

	print(paste("Finding Markers for ident :",i,":", sep=""))
	# find all markers of cluster i
	date()


	tryCatch({
		#	With lower thresholds?
		#	Occassionally, this errors
		#	Error in FindMarkers(object = ds, ident.1 = i, min.pct = 0.25) :
		#	  No genes pass logfc.threshold threshold
		#	Execution halted

		cluster.markers <- FindMarkers(object = ds, ident.1 = i, min.pct = 0.25)

		print(x = head(x = cluster.markers, n = 10))

		write.table(x = head(x = cluster.markers, n = 10),
			file=paste("genelist_cluster.",i,".csv",sep=""),
			col.names=NA,
			sep="\t")

	}, error = function(e) {
		#	Nothing. Just don't halt execution.
		print("Caught error.")
		print("Likely no genes meet the threshold.")
		print("Continuing.")
	})


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


print("Saving all environment")
save(list=ls(all=TRUE), file="Sample.RData")


print(paste0("Ending R script :",humanReadable(mem_used()),": at :", date(),":"))


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

