#!/usr/bin/env bash


set -e	#	exit if any command fails
set -u	#	Error on usage of unset variables
set -o pipefail



script=$( basename $0 )

#	Defaults:
#max=5
#num_cells=20000
genomedir="./myRef"
referencefasta="./myRef/myRef.fasta"

DROP_SEQ_PATH=~/.local/Drop-seq_tools-2.3.0
#PICARD_PATH=~/Downloads/picard.jar
STAR_PATH=~/.local/bin/STAR

function usage(){
	echo
	echo "Wrapper around calling Drop-seq_alignment.sh"
	echo
	echo "Script will loop over each bam file provided."
	echo
	echo "Usage:"
	echo
	echo "$script <OPTIONS> bam_file(s)"
	echo
	echo "Options:"
	echo "	--drop_seq /path/to/Drop-seq_alignment.sh : File or just directory"
	echo "	--star /path/to/STAR : File"
#	echo "	--picard /path/to/picard.jar : File or just directory"
#	echo "	--estimated-num-cells (-n) INTEGER : "
	echo "	--genomedir (-g) STRING : Directory of STAR genome directory"
	echo "	--referencefasta (-r) STRING : Reference fasta of the Drop-seq reference metadata bundle"
#	echo "	--max INTEGER ......... Maximum number of scans per file"
#	echo "	--verbose .................. NO VALUE, just boolean flag"
	echo
	echo "Default option values:"
	echo "	--drop_seq                ${DROP_SEQ_PATH}"
	echo "	--star                    ${STAR_PATH}"
#	echo "	--picard                  ${PICARD_PATH}"
#	echo "	--estimated-num-cells ... ${num_cells}"
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
#		-n|--estimated-num-cells)
#			shift; num_cells=$1; shift;;
		-g|--genomedir)
			shift; genomedir=$1; shift;;
		-r|--referencefasta)
			shift; referencefasta=$1; shift;;
		--drop_seq)
			shift; DROP_SEQ_PATH=$1; shift;;
		--star)
			shift; STAR_PATH=$1; shift;;
#		--picard)
#			shift; PICARD_PATH=$1; shift;;
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



if [ -f ${DROP_SEQ_PATH} ] ; then
	DROP_SEQ_PATH=$( dirname ${DROP_SEQ_PATH} )
fi





calling_dir=$PWD

while [ $# -ne 0 ] ; do
	cd $calling_dir
	echo "Processing :${1}:"

	bam_file_with_path=$1
#	bam_base=${bam_file_with_path%%.*}
#	bam_base=${bam_base##*/}
	bam_base=$( basename $1 .bam )
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

	if [ ! -d ${bam_base} ] ; then

		mkdir ${bam_base}

	fi

	echo "Running Drop Seq."

	if [ -f ${bam_base}/final.bam ] ; then
		echo "${bam_base}/final.bam exists. Drop seq already run."
	else

		if [ ! -f "${DROP_SEQ_PATH}/Drop-seq_alignment.sh" ] ; then
			echo "${DROP_SEQ_PATH}/Drop-seq_alignment.sh not found and is required."
			echo "Cannot continue."
			exit 2
		fi


#		if [ "${DROP_SEQ_VERSION}" == "2.0.0" ] ; then
			#	2.0.0 no longer accepts this option
			#			-n ${num_cells} \
			cmd="${DROP_SEQ_PATH}/Drop-seq_alignment.sh \
				-s ${STAR_PATH} \
				-g ${genomedir} \
				-r ${referencefasta} \
				-o ${bam_base} \
				${bam_file_with_path}"
			echo $cmd
			$cmd

#		else
#			#	Assuming version 1.13
#			cmd="${DROP_SEQ_PATH}/Drop-seq_alignment.sh \
#				-n ${num_cells} \
#				-g ${genomedir} \
#				-r ${referencefasta} \
#				-o ${bam_base} \
#				${bam_file_with_path}"
#			echo $cmd
#			$cmd

#			#	This is done simply to exercise the new naming conventions of 2.0.0
#			mv ${bam_base}/error_detected.bam ${bam_base}/final.bam
#			mv ${bam_base}/error_detected.bai ${bam_base}/final.bai

#		fi

	fi

	if [ ! -f "${bam_base}/final.bam" ] ; then
		echo "final.bam not found. Drop Seq failed? Cannot continue."
		exit 3
	fi


	cd "${bam_base}"
	echo $PWD
	ls -trail


	echo "Running BAMTagHistogram."

	if [ -f out_cell_readcounts.txt.gz ] ; then

		echo "out_cell_readcounts.txt.gz exists. BAMTagHistogram already run."

	else

#		if [ "${DROP_SEQ_VERSION}" == "2.0.0" ] ; then
			#	SERIOUSLY! BAM to Bam. I don't even think that Linux would let me link that. Rename during install?
			cmd="${DROP_SEQ_PATH}/BamTagHistogram \
				INPUT=final.bam \
				OUTPUT=out_cell_readcounts.txt.gz \
				TAG=XC"
#		else
#			#	Assuming version 1.13
#			cmd="${DROP_SEQ_PATH}/BAMTagHistogram \
#				INPUT=final.bam \
#				OUTPUT=out_cell_readcounts.txt.gz \
#				TAG=XC"
#		fi
		echo $cmd
		$cmd

	fi

	if [ -f out_cell_readcounts.txt.gz ] && [ -f cell_bc_file.txt.gz ] ; then

		echo "cell_bc_file.txt.gz already exists."

	else

		#
		#	Only keep those with more than one. This is just a test.
		#
		#zcat out_cell_readcounts.txt.gz | tail -n +2 | awk '( $1 > 1 ){print $2}' | gzip > cell_bc_file.txt.gz
		zcat out_cell_readcounts.txt.gz | tail -n +2 | awk '{print $2}' | gzip > cell_bc_file.txt.gz

	fi


#	2.3.0 producing ...
#	2.0.0 does something similar? I've run this before! How is it failing now!
#	/broad/hptmp/jake/sortingcollection.270241350814330195.tmp
#	Gotta set TMP_DIR


#********** The command line looks like this in the new syntax:
#**********
#**********    DigitalExpression -INPUT final.bam -OUTPUT final.dge.txt.gz -CELL_BC_FILE cell_bc_file.txt.gz -SUMMARY out_gene_exon_tagged.dge.summary.txt -MIN_NUM_GENES_PER_CELL 100
#**********
#
#
#15:14:23.478 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/home/jake/.local/Drop-seq_tools-2.3.0/jar/lib/picard-2.18.14.jar!/com/intel/gkl/native/libgkl_compression.so
#15:14:23.484 WARN  NativeLibraryLoader - Unable to load libgkl_compression.so from native/libgkl_compression.so (No such file or directory)
#15:14:23.485 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/home/jake/.local/Drop-seq_tools-2.3.0/jar/lib/picard-2.18.14.jar!/com/intel/gkl/native/libgkl_compression.so
#15:14:23.485 WARN  NativeLibraryLoader - Unable to load libgkl_compression.so from native/libgkl_compression.so (No such file or directory)
#[Fri Jul 26 15:14:23 PDT 2019] DigitalExpression SUMMARY=out_gene_exon_tagged.dge.summary.txt OUTPUT=final.dge.txt.gz INPUT=final.bam MIN_NUM_GENES_PER_CELL=100 CELL_BC_FILE=cell_bc_file.txt.gz    OUTPUT_READS_INSTEAD=false CELL_BARCODE_TAG=XC MOLECULAR_BARCODE_TAG=XM EDIT_DISTANCE=1 READ_MQ=10 MIN_BC_READ_THRESHOLD=0 USE_STRAND_INFO=true RARE_UMI_FILTER_THRESHOLD=0.0 GENE_NAME_TAG=gn GENE_STRAND_TAG=gs GENE_FUNCTION_TAG=gf STRAND_STRATEGY=SENSE LOCUS_FUNCTION_LIST=[CODING, UTR] VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
#[Fri Jul 26 15:14:23 PDT 2019] Executing as jake@system76-server on Linux 4.15.0-54-generic amd64; OpenJDK 64-Bit Server VM 1.8.0_212-8u212-b03-0ubuntu1.16.04.1-b03; Deflater: Jdk; Inflater: Jdk; Provider GCS is not available; Picard version: 2.3.0(34e6572_1555443285)
#INFO	2019-07-26 15:14:23	BarcodeListRetrieval	Found 54863 cell barcodes in file
#INFO	2019-07-26 15:14:23	DigitalExpression	Calculating digital expression for [54863] cells.
#15:14:23.639 WARN  IntelDeflaterFactory - IntelInflater is not supported, using Java.util.zip.Inflater
#[Fri Jul 26 15:14:37 PDT 2019] org.broadinstitute.dropseqrna.barnyard.DigitalExpression done. Elapsed time: 0.23 minutes.
#Runtime.totalMemory()=2990014464
#Exception in thread "main" htsjdk.samtools.util.RuntimeIOException: java.nio.file.NoSuchFileException: /jake/sortingcollection.679353111144896655.tmp
#	at htsjdk.samtools.util.SortingCollection.spillToDisk(SortingCollection.java:268)
#	at htsjdk.samtools.util.SortingCollection.add(SortingCollection.java:183)
#	at org.broadinstitute.dropseqrna.utils.SortingIteratorFactory.create(SortingIteratorFactory.java:67)
#	at org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory.create(SamRecordSortingIteratorFactory.java:57)
#	at org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator.<init>(UMIIterator.java:133)
#	at org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator.<init>(UMIIterator.java:70)
#	at org.broadinstitute.dropseqrna.barnyard.DigitalExpression.digitalExpression(DigitalExpression.java:188)
#	at org.broadinstitute.dropseqrna.barnyard.DigitalExpression.doWork(DigitalExpression.java:154)
#	at picard.cmdline.CommandLineProgram.instanceMain(CommandLineProgram.java:295)
#	at picard.cmdline.PicardCommandLine.instanceMain(PicardCommandLine.java:103)
#	at org.broadinstitute.dropseqrna.cmdline.DropSeqMain.main(DropSeqMain.java:42)
#Caused by: java.nio.file.NoSuchFileException: /jake/sortingcollection.679353111144896655.tmp
#	at sun.nio.fs.UnixException.translateToIOException(UnixException.java:86)
#	at sun.nio.fs.UnixException.rethrowAsIOException(UnixException.java:102)
#	at sun.nio.fs.UnixException.rethrowAsIOException(UnixException.java:107)
#	at sun.nio.fs.UnixFileSystemProvider.newByteChannel(UnixFileSystemProvider.java:214)
#	at java.nio.file.Files.newByteChannel(Files.java:361)
#	at java.nio.file.Files.createFile(Files.java:632)
#	at java.nio.file.TempFileHelper.create(TempFileHelper.java:138)
#	at java.nio.file.TempFileHelper.createTempFile(TempFileHelper.java:161)
#	at java.nio.file.Files.createTempFile(Files.java:852)
#	at htsjdk.samtools.util.IOUtil.newTempPath(IOUtil.java:328)
#	at htsjdk.samtools.util.SortingCollection.newTempFile(SortingCollection.java:279)
#	at htsjdk.samtools.util.SortingCollection.spillToDisk(SortingCollection.java:250)
#	... 10 more

	echo "Running DigitalExpression."

	if [ -f final.dge.txt.gz ] ; then

		echo "final.dge.txt.gz exists. DigitalExpression already run."

	else

		#	Found discussion on this. Said that it was fixed in 2.0.0, but apparently not.
		#	https://github.com/broadinstitute/Drop-seq/issues/61
		#	Passing a TMP_DIR fixes

		cmd="${DROP_SEQ_PATH}/DigitalExpression \
			TMP_DIR=/tmp/jake/ \
			INPUT=final.bam \
			OUTPUT=final.dge.txt.gz \
			CELL_BC_FILE=cell_bc_file.txt.gz \
			SUMMARY=out_gene_exon_tagged.dge.summary.txt \
			MIN_NUM_GENES_PER_CELL=100"
		echo $cmd
		$cmd

	fi


	if [ -f InitialSeuratObjectSample.RData ] ; then

		echo "InitialSeuratObjectSample.RData exists. create_seurat.R already run."

	else

		if [ -f final.dge.txt.gz ] ; then

			#	If its too big, it will run for days and then run out of memory/swap space
			#70514494 B6/final.dge.txt.gz
			if [ $( wc -c < final.dge.txt.gz ) -lt 70000000 ] ; then

				echo "Running create_seurat.R"

				cmd=create_seurat.R
				echo $cmd
				$cmd

			fi

		fi

	fi


	if [ -f InitialSeuratObjectSample.RData ] ; then

		if [ -f Rplots.pdf ] ; then

			echo "Rplots.pdf exists. seurat.R already run."

		else

			echo "Running seurat.R"

			#	R is pretty bad at garbage collection.
			#	Reading final.dge.txt.gz and creating the seurat object then quiting.
			#	Then running another script that reads in the seurat works well.

			#		cmd="seurat.R --redo"	#	this is the default now
			cmd="seurat.R"
			echo $cmd
			$cmd

		fi

	fi

	echo
	shift
done

