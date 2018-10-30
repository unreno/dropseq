#!/usr/bin/env bash


script=`basename $0`

#	Defaults:
#max=5
num_cells=20000
genomedir="$HOME/working/mm10c_star"
referencefasta="$HOME/mm10c/mm10c.fasta"

function usage(){
	echo
	echo "Wrapper around calling Drop-seq_alignment.sh"
	echo
	echo "Usage:"
	echo
	echo "$script <OPTIONS> bam_file(s)"
	echo
	echo "Options:"
	echo "	--estimated-num-cells (-n) INTEGER : "
	echo "	--genomedir (-g) STRING : Directory of STAR genome directory"
	echo "	--referencefasta (-r) STRING : Reference fasta of the Drop-seq reference metadata bundle"
#	echo "	--max INTEGER ......... Maximum number of scans per file"
#	echo "	--verbose .................. NO VALUE, just boolean flag"
	echo
	echo "Default option values:"
	echo "	--estimated-num-cells ... ${num_cells}"
	echo "	--genomedir ............. ${genomedir}"
	echo "	--referencefasta ........ ${referencefasta}"
#	echo "	--max ........ ${max}"
#	echo "	--verbose .......... ${verbose}"
	echo
	echo "Examples:"
	echo "	$script file.bam"
#	echo "	$script --max 10 file.mgf"
	echo
	echo
	exit 1
}

while [ $# -ne 0 ] ; do
	case $1 in
		-n|--estimated-num-cells)
			shift; num_cells=$1; shift;;
		-g|--genomedir)
			shift; genomedir=$1; shift;;
		-r|--referencefasta)
			shift; referencefasta=$1; shift;;
#		-m|--m*)
#			shift; max=$1; shift;;
#		-v|--v*)
#			verbose=true; shift;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done

#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage






calling_dir=$PWD

while [ $# -ne 0 ] ; do
	cd $calling_dir
	echo "Processing :${1}:"

	bam_file_with_path=$1
	bam_base=${bam_file_with_path%%.*}
	bam_base=${bam_base##*/}
	mkdir -p "${bam_base}"
	
	#	prototype script for AWS AMI so many hard coded values
	#		
	#		-g <genomedir>      : Directory of STAR genome directory.  Required.
	#		-r <referencefasta> : Reference fasta of the Drop-seq reference metadata bundle.  Required.
	#		-n <estimated-num-cells> : estimate of number of cells in experiment.  Required.
	#		-d <dropseq_root>   : Directory containing Drop-seq executables.  Default: directory containing this script.
	#		-o <outputdir>      : Where to write output bam.  Default: current directory.
	#		-t <tmpdir>         : Where to write temporary files.  Default: a new subdirectory in .
	#		-s <STAR_path>      : Full path of STAR.  Default: STAR is found via PATH environment variable.
	#		-p                  : Reduce file I/O by pipeline commands together.  Requires more memory and processing power.
	#		-e                  : Echo commands instead of executing them.  Cannot use with -p.
	#		

	#	Last line of Drop-seq_alignment deletes all tmp files, so I commented that line out to see if useful.
	#	Should be an option

	cmd="Drop-seq_alignment.sh \
		-g ${genomedir} \
		-r ${referencefasta} \
		-n ${num_cells} \
		-o ${bam_base} \
		${bam_file_with_path}"
	echo $cmd
	$cmd

	cd "${bam_base}"


	BAMTagHistogram \
		INPUT=error_detected.bam \
		OUTPUT=out_cell_readcounts.txt.gz \
		TAG=XC

	#
	#	Only keep those with more than one. This is just a test.
	#
	#zcat out_cell_readcounts.txt.gz | tail -n +2 | awk '( $1 > 1 ){print $2}' | gzip > cell_bc_file.txt.gz
	zcat out_cell_readcounts.txt.gz | tail -n +2 | awk '{print $2}' | gzip > cell_bc_file.txt.gz

	DigitalExpression \
		INPUT=error_detected.bam \
		OUTPUT=error_detected.dge.txt.gz \
		CELL_BC_FILE=cell_bc_file.txt.gz \
		SUMMARY=out_gene_exon_tagged.dge.summary.txt \
		MIN_NUM_GENES_PER_CELL=100

	create_seurat.R

	#	R is pretty bad at garbage collection.
	#	Reading error_detected.dge.txt.gz and creating the seurat object then quiting.
	#	Then running another script that reads in the seurat works well.

	seurat.R --redo

	echo
	shift
done

