# Drop Seq Processing

This repository is a collection of scripts for using and analysing Drop Seq data.

It was developed with older version but currently runs with the following requirements.

* STAR - STAR-2.7.1a.tar.gz ( https://github.com/alexdobin/STAR/releases )
  * wget -O STAR-2.7.1a.tar.gz https://github.com/alexdobin/STAR/archive/2.7.1a.tar.gz
* Drop Seq 2.3.0
  * Drop\_seq - Drop-seq\_tools-2.3.0.zip ( http://mccarrolllab.com/dropseq/ )
  * wget https://github.com/broadinstitute/Drop-seq/releases/download/v2.3.0/Drop-seq_tools-2.3.0.zip
* R version 3.5.1 (2018-07-02) -- "Feather Spray"
* R packages
  * if (!requireNamespace("BiocManager", quietly = TRUE))
    * install.packages("BiocManager")
  * BiocManager::install(c("devtools","Seurat","pryr","gdata","dplyr"),update = TRUE, ask = FALSE)

* R packages ( Seurat 3 is much different than Seurat 2. )
  * require(devtools)
  * install_version('Seurat',version='2.3.4')






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

Drop Seq requires a reference set with matching `.fasta`, `.gtf`, `.dict`, `.refFlat` and `STAR` database.

I believe that the 4 reference files must have the same root name. (e.q. mm10.fasta, mm10.gtf, etc.)

For this example tutorial, we will be using a modified mm10 reference.


###	Create, or obtain, a required reference `.fasta` file


We can start by downloading the complete mm10 reference metadata from NIH.

ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63472/suppl/GSE63472_mm10_reference_metadata.tar.gz

This `.tar.gz` contains all of the necessary files, except the actual STAR database itself.


If that is not available, you could also download a compilation of this reference in many formats.

ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz

https://support.illumina.com/sequencing/sequencing_software/igenome.html

It contains the raw `.fasta`, `.dict`, `.gtf`, and `.refFlat` which would likely need renamed to "mm10" before use.

It also contains prepared indexes for aligners `bwa`, `bowtie` and `bowtie2` as well as other products.


You could also get the individual `.fasta.gz` files from http://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/ and combine them. Of course, this doesn't include the `.gtf` file.


If you have any additional sequences you'd like to add to your reference, now is the time.
Append this existing `.fasta` file with these sequences for your reference.

```
>eGFP
ATGGTGAGCAAGGGC ...
>SV40polya
ATCTAGATAACTGAT ...
```

###	Create required reference `.gtf` file

We'll be using and appending the `mm10.gtf` file.

If you have any additional sequences you'd like to add to your reference, you will also need to add them to your `.gtf` file.
Append this existing `.gtf` file for use by `STAR` like so ...

```
eGFP	AddedGenes	exon	1	720	.	+	0	gene_id "eGFP"; gene_name "eGFP"; transcript_id "eGFP"; transcript_name "eGFP";
SV40polya	AddedGenes	exon	1	240	.	+	0	gene_id "SV40polya"; gene_name "SV40polya"; transcript_id "SV40polya"; transcript_name "SV40polya";
```

A `.gtf` file is a TAB separated file so edit with caution.

JUST for your viewing pleasure, the TABs have been converted to PIPEs. REMINDER, a `.gtf` file is a TAB separated file. Don't use PIPEs in the actual file.

```
eGFP|AddedGenes|exon|1|720|.|+|0|gene_id "eGFP"; gene_name "eGFP"; transcript_id "eGFP"; transcript_name "eGFP";
SV40polya|AddedGenes|exon|1|240|.|+|0|gene_id "SV40polya"; gene_name "SV40polya"; transcript_id "SV40polya"; transcript_name "SV40polya";
```

That final column MUST contain the subfields `gene_id`, `gene_name`, `transcript_id`, and `transcript_name`.








###	Compare / Check `.fasta` and `.gtf`


**The `.fasta` and `.gtf` sequence names MUST, MUST, MUST be consistent.**

`GSE63472_mm10_reference_metadata.tar.gz` IS NOT!

The `.fasta` file includes the "chr" prefix, whilst the `.gtf` does not. Why?

I editted the `.fasta` stripping the "chr" prefix before recreating the `.dict` and `.refFlat` files.

( 20190726 - I just checked this and it is fine. Neither have "chr". Must've had a different copy.)


This isn't the best check, but its something.
There needs to be some overlap in sequence names between the 2 files.
In this case there are 24.

```BASH
grep "^>" mm10c.fasta | sed 's/>//' | sort > mm10c.fasta.sequences
awk -F"\t" '{print $1}' mm10c.gtf | sort | uniq > mm10c.gtf.sequences
comm -12 mm10c.fasta.sequences mm10c.gtf.sequences | wc -l
24
```





###	Create required reference `.dict` from `.fasta` using Picard's CreateSequenceDictionary

This only takes a few seconds to generate.

If you haven't modified the reference `.fasta`, this is unnecessary.

```BASH
export DROP_SEQ_PATH=~/Drop-seq_tools-1.13
export PICARD_PATH=${DROP_SEQ_PATH}/3rdParty/picard/
java -jar $PICARD_PATH/picard.jar CreateSequenceDictionary REFERENCE=mm10/mm10.fasta
```

The `.dict` file should have 1 more line that the number of sequences in the reference `.fasta` file.
It begins with an @HD line and is followed by an @SQ line for each sequence.
It is very similar to a `.sam` file header. 

For example, 

```BASH
wc -l mm10.dict 
69 mm10.dict

grep -c "^@SQ" mm10.dict 
68

grep -c "^>" mm10.fasta 
68
```


###	Create required reference `.refFlat` file

If you haven't modified the reference `.fasta` and you already have a `.dict` file, this is unnecessary.
Not too long to create here either.

Using Drop Seq tools' `ConvertToRefFlat`

```BASH
export DROP_SEQ_PATH=~/Drop-seq_tools-1.13
$DROP_SEQ_PATH/ConvertToRefFlat ANNOTATIONS_FILE=mm10/mm10.gtf SEQUENCE_DICTIONARY=mm10/mm10.dict OUTPUT=mm10/mm10.refFlat
```



I get a lot of "GTFReader	Multiple gene IDs for gene".

This is odd as it is using the `GSE63472_mm10_reference_metadata.tar.gz` reference `.gtf` and `.dict` when trying to recreate the provided `.refFlat`.



Alternatively, the utility `gtfToGenePred` also creates a `.refFlat` file.
I obtained this program from UCSC's Kent Utils.
It does require some options and does produce a `.refFlat` with columns in a different order so some post run manipulation is required as well.
I do find it a bit unacceptable that this gene data format is not standardized.


```BASH
gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons mm10.gtf /dev/stdout | awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'
```

Of course, you can also manually edit this file to include any of your additions.

```BASH
SV40polya	SV40polya	SV40polya	+	0	240	240	240	1	0,	240,
eGFP	eGFP	eGFP	+	0	720	720	720	1	0,	720,
```


The `.refFlat` file should have a line for each of the transcript ids in the `.gtf` file.
Minus any removed due to error checks.


For example,

```BASH
wc -l mm10.gtf mm10.refFlat 
  1131185 mm10.gtf
    88228 mm10.refFlat
  1219413 total

sed -e 's/^.*transcript_id "//' -e 's/".*$//' mm10.gtf > mm10.gtf.transcript_ids

sort mm10.gtf.transcript_ids | uniq | wc -l
88228
```


###	Create actual STAR reference

And finally, create a `STAR` reference database.
This can take a couple hours.
Create a subdir for 8 or so files generated here.
There are other options you can explore.

```BASH
mkdir mm10STAR
STAR --genomeFastaFiles mm10/mm10.fasta --runMode genomeGenerate --genomeDir $PWD/mm10STAR --sjdbGTFfile mm10/mm10.gtf --runThreadN 40
```

Note:

If the reference fasta is small, or if you are getting seg faults during the alignment, you may wish to recreate the index with a lower the value for option `--genomeSAindexNbases`. Its default is 14, so try 5 or something. On a previous run, my index only had 2 sequences in it and this fixed it for me.






##	Use a Docker Image instead

The following will create the above as a docker image with STAR 2.5.3a, Drop Seq tools 1.13 and a modified mm10 reference.

Be advised, building this image will take a couple hours and require a machine with at least 20GB of memory to make the STAR reference.


```BASH
mkdir tmp; cd tmp
wget https://github.com/unreno/dropseq/blob/master/docker/Dockerfile-1.13
docker build -t dropseq1 -f Dockerfile-1.13 .
docker run -it dropseq1
```

Of course, you can do this whereever and call the image whatever you like.

An docker file to create an image with STAR 2.6.1c and Drop Seq tools 2.0.0 is also available.

```BASH
mkdir tmp; cd tmp
wget https://github.com/unreno/dropseq/blob/master/docker/Dockerfile-2.0.0
docker build -t dropseq2 -f Dockerfile-2.0.0 .
docker run -it dropseq2
```







##	Prepare your dataset

###	Combine fastq files into unaligned bam file

Drop Seq's script expects a sorted, unaligned bam as primary input.

You can do this with Picard's FastqToSam tool, if you have fastq files.

```BASH
export DROP_SEQ_PATH=~/Drop-seq_tools-1.13
export PICARD_PATH=${DROP_SEQ_PATH}/3rdParty/picard/
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


export DROP_SEQ_PATH=~/Drop-seq_tools-1.13
export PICARD_PATH=${DROP_SEQ_PATH}/3rdParty/picard/
java -jar $PICARD_PATH/picard.jar FastqToSam \
	F1=B4_S2_L001_R1_001.fastq.gz \
	F2=B4_S2_L001_R2_001.fastq.gz \
	O=B4_S2_L001.bam SM=B4_S2_L001
java -jar $PICARD_PATH/picard.jar FastqToSam \
	F1=B4_S2_L002_R1_001.fastq.gz \
	F2=B4_S2_L002_R2_001.fastq.gz \
	O=B4_S2_L002.bam SM=B4_S2_L002
java -jar $PICARD_PATH/picard.jar FastqToSam \
	F1=B4_S2_L003_R1_001.fastq.gz \
	F2=B4_S2_L003_R2_001.fastq.gz \
	O=B4_S2_L003.bam SM=B4_S2_L003
java -jar $PICARD_PATH/picard.jar FastqToSam \
	F1=B4_S2_L004_R1_001.fastq.gz \
	F2=B4_S2_L004_R2_001.fastq.gz \
	O=B4_S2_L004.bam SM=B4_S2_L004

```


###	Merge sample bam files

If your sample is comprised of multiple pairs of FASTQ files, merge them with something like ...

```BASH
export DROP_SEQ_PATH=~/Drop-seq_tools-1.13
export PICARD_PATH=${DROP_SEQ_PATH}/3rdParty/picard/
java -jar $PICARD_PATH/picard.jar MergeSamFiles \
	INPUT=B3_S1_L001.bam \
	INPUT=B3_S1_L002.bam \
	INPUT=B3_S1_L003.bam \
	INPUT=B3_S1_L004.bam \
	ASSUME_SORTED=true \
	SORT_ORDER=queryname \
	OUTPUT=B3.bam


export DROP_SEQ_PATH=~/Drop-seq_tools-1.13
export PICARD_PATH=${DROP_SEQ_PATH}/3rdParty/picard/
java -jar $PICARD_PATH/picard.jar MergeSamFiles \
	INPUT=B4_S2_L001.bam \
	INPUT=B4_S2_L002.bam \
	INPUT=B4_S2_L003.bam \
	INPUT=B4_S2_L004.bam \
	ASSUME_SORTED=true \
	SORT_ORDER=queryname \
	OUTPUT=B4.bam
```



##	Run Drop Seq Wrapper


Be advised, portions of this pipeline require a GREAT deal of memory.
If it is unavailable, the script will simply be killed and crash.


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
	--drop_seq                /Users/jakewendt/Drop-seq_tools-1.13
	--estimated-num-cells ... 20000
	--genomedir ............. ./myRef
	--referencefasta ........ ./myRef/myRef.fasta

Examples:
	drop_seq.bash file.bam

```


Drop Seq's script expects an unaligned bam as primary input.



```BASH
export DROP_SEQ_PATH=~/Drop-seq_tools-1.13

drop_seq.bash \
	--drop_seq ${DROP_SEQ_PATH}/Drop-seq_alignment.sh \
	--estimated-num-cells 20000 \
	--genomedir ${PWD}/myRefSTAR \
	--referencefasta ${PWD}/myRef/myRef.fasta \
	B3.bam
```

```BASH
export DROP_SEQ_PATH=~/Drop-seq_tools-1.13

nohup drop_seq.bash --drop_seq ${DROP_SEQ_PATH}/Drop-seq_alignment.sh \
	--estimated-num-cells 20000 --genomedir ${PWD}/mm10STAR \
	--referencefasta ${PWD}/mm10/mm10.fasta B3_S1_L001.bam > B3_S1_L001.log &

nohup drop_seq.bash --drop_seq ${DROP_SEQ_PATH}/Drop-seq_alignment.sh \
	--estimated-num-cells 20000 --genomedir ${PWD}/mm10STAR \
	--referencefasta ${PWD}/mm10/mm10.fasta B3.bam > B3.log &

nohup drop_seq.bash --drop_seq ${DROP_SEQ_PATH}/Drop-seq_alignment.sh \
	--estimated-num-cells 20000 --genomedir ${PWD}/mm10STAR \
	--referencefasta ${PWD}/mm10/mm10.fasta B4.bam > B4.log &
```


##	Run Drop Seq on docker image


Create a directory to share with the docker image, place you sample bam in it and execute the docker image.

```BASH
mkdir ~/local_to_share
cp SAMPLE.bam ~/local_to_share

docker run -v ~/local_to_share:/shared --rm dropseq1 bash -c "cd /shared/; drop_seq.bash --drop_seq /root/Drop-seq_tools --genomedir /root/mm10STAR --referencefasta /root/mm10/mm10.fasta SAMPLE.bam > SAMPLE.log 2>&1"

tail --retry -f ~/local_to_share/SAMPLE.log
```

Of course, you could just start a docker container and be more interactive.

```BASH
docker run -it dropseq1
```

