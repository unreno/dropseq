#!/usr/bin/env bash



set -e	#	exit if any command fails
set -u	#	Error on usage of unset variables
set -o pipefail


export DROP_SEQ_PATH=~/.local/Drop-seq_tools-2.3.0
export PICARD_PATH=${DROP_SEQ_PATH}/3rdParty/picard/


f=GSE63472_mm10_reference_metadata.tar.gz
if [ -f $f ] && [ ! -w $f ] ; then
	echo "Write-protected $f exists. Skipping."
else
	echo "Creating $f"
	rm ${f}
	wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63472/suppl/$f
	chmod a-w $f
fi

f=mm10
if [ -d $f ] && [ ! -w $f ] ; then
	echo "Write-protected $f exists. Skipping."		#	DIRECTORY
else
	echo "Creating $f"
	tar xfvz GSE63472_mm10_reference_metadata.tar.gz
	chmod a-x $f/*
	chmod -R a-w $f
fi


mkdir -p mm10x

for ext in fasta gtf dict refFlat ; do

	f=mm10x/mm10x.${ext}
	if [ -f $f ] && [ ! -w $f ] ; then
		echo "Write-protected $f exists. Skipping."
	else
		echo "Creating $f"
		cat mm10/mm10.${ext} ../misc/*${ext} > ${f}
		chmod a-w $f
	fi

done




#	I'm not sure if this is a good idea or not, but I have no idea of the expected
#	location of the 2 added genes so just copying the original mm10 files.
cp mm10/mm10.exons.intervals mm10x/mm10x.exons.intervals
cp mm10/mm10.genes.intervals mm10x/mm10x.genes.intervals







#f=mm10xSTAR
#mkdir -p ${f}
#if [ -d $f ] && [ ! -w $f ] ; then
#	echo "Write-protected $f exists. Skipping."		#	DIRECTORY
#else
#	echo "Creating $f"
#	~/.local/STAR-2.7.1a/bin/Linux_x86_64/STAR \
#		--genomeFastaFiles mm10x/mm10x.fasta \
#		--runMode genomeGenerate \
#		--genomeDir $PWD/mm10xSTAR \
#		--sjdbGTFfile mm10x/mm10x.gtf \
#		--runThreadN 40
#	chmod -R a-w $f
#fi

#	When running drop seq on 2.7.1c, I get ...
#	EXITING because of FATAL ERROR: Genome version is INCOMPATIBLE with current STAR version
#	SOLUTION: please re-generate genome from scratch with the latest version of STAR
#	which is seems odd as I just did this?

#	My bad. I let DropSeq use an old version of STAR in my PATH.


for star in 2.7.1a 2.6.1c 2.5.3a ; do
	f=mm10x.STAR-${star}
	mkdir -p ${f}
	if [ -d $f ] && [ ! -w $f ] ; then
		echo "Write-protected $f exists. Skipping."		#	DIRECTORY
	else
		echo "Creating $f"
		~/.local/STAR-${star}/bin/Linux_x86_64/STAR \
			--genomeFastaFiles mm10x/mm10x.fasta \
			--runMode genomeGenerate \
			--genomeDir $PWD/${f} \
			--sjdbGTFfile mm10x/mm10x.gtf \
			--runThreadN 40
		chmod -R a-w $f
	fi
done

