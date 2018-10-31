# Drop Seq Processing

This repository is a collection of scripts for using and analysing Drop Seq data.

It was developed with the following requirements.

* STAR - STAR-2.5.3a.tar.gz  ( https://github.com/alexdobin/STAR/releases )
* Picard - picard.jar ( http://broadinstitute.github.io/picard/ )
* Drop\_seq - Drop-seq\_tools-1.13-3.zip ( http://mccarrolllab.com/dropseq/ )
* R version 3.5.1 (2018-07-02) -- "Feather Spray"
  * source("https://bioconductor.org/biocLite.R")
  * biocLite( c("devtools","Seurat","pryr","gdata","optparse") )


Newer versions may work, but they have not been tested.


##	Installation

* Download this repository.
* Either ...
  * Copy or link the Makefile.example to Makefile
  * Edit the Makefile as necessary to install scripts where you'd like
  * Run `make install`
* ... or ...
  * Copy the scripts from scripts/ to where ever you want them
* ... or ...
  * Run the scripts in scripts/ from where they are




##	Prepare STAR reference


###	Create reference ".fasta" file

Perhaps download and combine the fasta.gz files from http://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/

Append existing or create a new fasta file with any special sequences for your reference.

```
>eGFP
ATGGTGAGCAAGGGC ...
>SV40polya
ATCTAGATAACTGAT ...
```

###	Create reference ".gtf" file

Append existing or create a new .gtf file for use by STAR like so ...

```
eGFP	AddedGenes	exon	1	720	.	+	0	gene_id "eGFP"; gene_name "eGFP"; transcript_id "eGFP"; transcript_name "eGFP";
SV40polya	AddedGenes	exon	1	240	.	+	0	gene_id "SV40polya"; gene_name "SV40polya"; transcript_id "SV40polya"; transcript_name "SV40polya";
```

Or like so where the tabs have been converted to pipes (JUST for your viewing pleasure. A gtf file is a TAB separated file.)

```
eGFP|AddedGenes|exon|1|720|.|+|0|gene_id "eGFP"; gene_name "eGFP"; transcript_id "eGFP"; transcript_name "eGFP";
SV40polya|AddedGenes|exon|1|240|.|+|0|gene_id "SV40polya"; gene_name "SV40polya"; transcript_id "SV40polya"; transcript_name "SV40polya";
```


###	And create a reference ".dict" using Picard's CreateSequenceDictionary

This only takes a few seconds

```BASH
export PICARD_PATH=~/Downloads/
java -jar $PICARD_PATH/picard.jar CreateSequenceDictionary REFERENCE=myRef/myRef.fasta
```

###	Create reference ".refFlat" file

Not too long here either.

Using Drop Seq tools' ConvertToRefFlat

```BASH
export DROP_SEQ_PATH=~/Downloads/Drop-seq_tools-1.13
$DROP_SEQ_PATH/ConvertToRefFlat ANNOTATIONS_FILE=myRef/myRef.gtf SEQUENCE_DICTIONARY=myRef/myRef.dict OUTPUT=myRef/myRef.refFlat
```

###	Create actual STAR reference

```BASH
mkdir myRefSTAR
STAR --genomeFastaFiles myRef/myRef.fasta --runMode genomeGenerate --genomeDir $PWD/myRefSTAR --sjdbGTFfile myRef/myRef.gtf 
```

If the reference fasta is small, or if you are getting seg faults during the alignment, you may wish to recreate the index with a lower the value for option `--genomeSAindexNbases`. Its default is 14, so try 5 or something. My index only had 2 sequences in it and this fixed it for me.








##	Run Drop Seq Wrapper

###	Prepare your dataset

####	Combine fastq files into unaligned bam file

Drop Seq's script expects an unaligned bam as primary input.

```BASH
java -jar $PICARD_PATH/picard.jar FastqToSam \
	F1=B3_S1_L001_R1_001.fastq.gz \
	F2=B3_S1_L001_R2_001.fastq.gz \
	O=B3_S1_L001.bam SM=B3_S1_L001

java -jar $PICARD_PATH/picard.jar FastqToSam \
	F1=B3_S1_L002_R1_001.fastq.gz \
	F2=B3_S1_L002_R2_001.fastq.gz \
	O=B3_S1_L002.bam SM=B3_S1_L002

java -jar $PICARD_PATH/picard.jar FastqToSam \
	F1=B3_S1_L003_R1_001.fastq.gz \
	F2=B3_S1_L003_R2_001.fastq.gz \
	O=B3_S1_L003.bam SM=B3_S1_L003

java -jar $PICARD_PATH/picard.jar FastqToSam \
	F1=B3_S1_L004_R1_001.fastq.gz \
	F2=B3_S1_L004_R2_001.fastq.gz \
	O=B3_S1_L004.bam SM=B3_S1_L004
```


####	Merge sample bam files

If your sample is comprised of multiple pairs of FASTQ files, merge them with ...

```BASH
java -jar $PICARD_PATH/picard.jar MergeSamFiles \
	INPUT=B3_S1_L001.bam \
	INPUT=B3_S1_L002.bam \
	INPUT=B3_S1_L003.bam \
	INPUT=B3_S1_L004.bam \
	ASSUME_SORTED=true \
	SORT_ORDER=queryname \
	OUTPUT=B3.bam
```




##	Execute

```BASH
drop_seq.bash

Wrapper around calling Drop-seq_alignment.sh

Script will loop over each bam file provided.

Usage:

drop_seq.bash <OPTIONS> bam_file(s)

Options:
	--drop_seq /path/to/Drop-seq_alignment.sh : File or just directory
	--estimated-num-cells (-n) INTEGER : 
	--genomedir (-g) STRING : Directory of STAR genome directory
	--referencefasta (-r) STRING : Reference fasta of the Drop-seq reference metadata bundle

Default option values:
	--drop_seq                /Users/jakewendt/Downloads/Drop-seq_tools-1.13
	--estimated-num-cells ... 20000
	--genomedir ............. ./myRef
	--referencefasta ........ ./myRef/myRef.fasta

Examples:
	drop_seq.bash file.bam

```


Drop Seq's script expects an unaligned bam as primary input.



```BASH
drop_seq.bash \
	--drop_seq ${DROP_SEQ_PATH}/Drop-seq_alignment.sh \
	--estimated-num-cells 20000 \
	--genomedir ${PWD}/myRefSTAR \
	--referencefasta ${PWD}/myRef/myRef.fasta \
	B3.bam
```




##	Analysis

If your sample or your reference is too small, the output will be small and the following analysis will likely fail due to zero's in the data.


```BASH
R
```

Run in two steps because R is bad at memory management.

```BASH
create_seurat.R

seurat.R --redo
```

