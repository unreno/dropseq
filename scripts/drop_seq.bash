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
	
#	tmp=${bam_base}/tmp
#	mkdir -p "${tmp}"

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

#	Once STAR in path, can remove from here.
#	Once DropSeq in path, can remove path from command (NEED TO MOVE (or link) ALL DROP SEQ SCRIPTS TO SAME DIR)
#	Add options for mm10 star ref dir and mm10 fasta dir

#	cmd="~/Drop-seq_tools-1.13/Drop-seq_alignment.sh \
#		-s ~/STAR-2.5.3a/bin/Linux_x86_64/STAR \
#		-g ~/working/mm10_star/ \
#		-r ~/mm10/mm10.fasta \
#		-t "${tmp}" \


#	EC2 no longer happy with ~ and absence of absolute STAR dir
#	may have to re-add this, and with /home/ec2-user/ instead of ~/
#
#	Seems unneeded now
#		-s ~/STAR-2.5.3a/bin/Linux_x86_64/STAR \

	cmd="Drop-seq_alignment.sh \
		-g ${genomedir} \
		-r ${referencefasta} \
		-n ${num_cells} \
		-o ${bam_base} \
		${bam_file_with_path}"
	echo $cmd
	$cmd

	cd "${bam_base}"
	cmd="dge.bash"
	echo $cmd
	$cmd

	echo
	shift
done
