#!/usr/bin/env bash



#	find . -name error_detected.bam -execdir ~/syryu/drop_seq_alignment/dge.bash \; > dge.log 2>&1


echo
echo
echo
pwd
echo



#~/singlecell/Drop-seq_tools-1.13/BAMTagHistogram \
BAMTagHistogram \
	INPUT=error_detected.bam \
	OUTPUT=out_cell_readcounts.txt.gz \
	TAG=XC

#
#	Only keep those with more than one. This is just a test.
#
#zcat out_cell_readcounts.txt.gz | tail -n +2 | awk '( $1 > 1 ){print $2}' | gzip > cell_bc_file.txt.gz
zcat out_cell_readcounts.txt.gz | tail -n +2 | awk '{print $2}' | gzip > cell_bc_file.txt.gz

#~/singlecell/Drop-seq_tools-1.13/DigitalExpression \
DigitalExpression \
	INPUT=error_detected.bam \
	OUTPUT=error_detected.dge.txt.gz \
	CELL_BC_FILE=cell_bc_file.txt.gz \
	SUMMARY=out_gene_exon_tagged.dge.summary.txt \
	MIN_NUM_GENES_PER_CELL=100

#	20180518 - set to 100 for latest data
#	MIN_NUM_GENES_PER_CELL=100

#	20180314 - refined for B3 and B4, but also to be reused for 1 and 2
#	MIN_NUM_GENES_PER_CELL=100

#	20180302 
#	undo NUM_CORE_BARCODES=100 MIN_NUM_GENES_PER_CELL=200 

#	20180301
#	was MIN_NUM_GENES_PER_CELL=100

#	20180122
#	was MIN_NUM_GENES_PER_CELL=100 NUM_CORE_BARCODES=100


#	NUM_CORE_BARCODES=100 ... Doesn't seem to make any difference. (Maybe default value?)

create_seurat.R

#	R is pretty bad at garbage collection.
#	Reading error_detected.dge.txt.gz and creating the seurat object then quiting.
#	Then running another script that reads in the seurat works well.

seurat.R --redo

echo
