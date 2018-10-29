
#	DropSeq Processing

##	20171227


###	Downloaded

* STAR - STAR-2.5.3a.tar.gz  ( https://github.com/alexdobin/STAR/releases )
* Picard - picard.jar ( http://broadinstitute.github.io/picard/ )
* Drop\_seq - Drop-seq\_tools-1.13-3.zip ( http://mccarrolllab.com/dropseq/ )
* Mouse reference files ... GSE63472\_mm10\_reference\_metadata.tar.gz


###	Converted all fastq files to unaligned bams with picard ( my script convert\_fastq\_files\_to\_bams.bash )

```
java -jar $basedir/picard.jar FastqToSam \
	F1=$fastq1 \
	F2=$fastq2 \
	O=$basename.bam \
	SM=$basename
```


###	Generated STAR reference from mm10.fasta

```
mkdir mm10_star
~/singlecell/STAR-2.5.3a/bin/MacOSX_x86_64/STAR --genomeFastaFiles mm10.fasta  --runMode genomeGenerate --genomeDir $PWD/mm10_star --runThreadN 8 --sjdbGTFfile mm10.gtf
```

This seemed to do nothing for a while, but used quite a bit of the CPU.
Then started producing files.
Eventually it crashed on my laptop.

Update Amazon AMI to include needed software and reference files.


```
create_ec2_instance.bash --profile syryu --key ~/.aws/JakeSYRyu.pem --instance-type c4.large 
sudo yum update
wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz
wget https://github.com/broadinstitute/picard/releases/download/2.17.0/picard.jar
wget -O Drop-seq_tools.zip http://mccarrolllab.com/download/1276/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63472/suppl/GSE63472_mm10_reference_metadata.tar.gz
```

removed java-1.7.0 and installed java-1.8.0 (picard needed)

created new ami with 12GB volume as needed more space.

c4.large (2,8,3.75)
  what():  std::bad\_alloc


Need more memory. Trying m4.xlarge (4,13,16) with --volume-size 100
4 vCPU, 13 ECU, 16GB memory. Trying with 13 threads. (shoulda used 4 instead)
It seems that STAR will "continue" if the files are there. AWESOME SAUCE!

```
mkdir ~/working/mm10_star
cd ~/working/
~/STAR-2.5.3a/bin/Linux_x86_64/STAR --genomeFastaFiles ~/mm10/mm10.fasta  --runMode genomeGenerate --genomeDir ~/working/mm10_star --runThreadN 13 --sjdbGTFfile ~/mm10/mm10.gtf
```

Installed htop. STAR PEGS the processors. Runs out of memory. Then crashes! (even with just 2 threads)
Trying m4.2xlarge (8,26,32)

```
create_ec2_instance.bash --profile syryu --key ~/.aws/JakeSYRyu.pem --instance-type m4.2xlarge --volume-size 100 --NOT-DRY-RUN

cd ~/working/
~/STAR-2.5.3a/bin/Linux_x86_64/STAR --genomeFastaFiles ~/mm10/mm10.fasta  --runMode genomeGenerate --genomeDir ~/working/mm10_star --runThreadN 8 --sjdbGTFfile ~/mm10/mm10.gtf
```

Pegged 8 processors. Memory seems to have stopped at 23GB. Phew. 26GB! Ehh. 27GB!! Up and down. Still going.
Done in about an hour. Not too bad. Final reference is about 25GB so rather than store, I'll regenerate if needed.



###	Copy up bam file for test run.


`scp -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no  /Users/jakewendt/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_03_L003-ds.ee0be2f1073b47d083679a1b81232597/2A_S2_L003.bam ec2-user@$ip:working/`




###	Running drop-seq's alignment script ...

```
mkdir ~/working/2A_S2_L003.dropseq.4000
cd ~/working/2A_S2_L003.dropseq.4000
~/Drop-seq_tools-1.13/Drop-seq_alignment.sh -g ~/working/mm10_star/ -s ~/STAR-2.5.3a/bin/Linux_x86_64/STAR -n 4000 -r ~/mm10/mm10.fasta ~/working/2A_S2_L003.bam
tar cfv - 2A_S2_L003.dropseq.4000 | gzip --best > 2A_S2_L003.dropseq.4000.tar.gz
```

GUESSED
-n <estimated-num-cells> : estimate of number of cells in experiment.  Required.  ( ~4000 from researcher )



This generated a lot of error bam data. (Apparently this is good.)


So Young recommends testing with ... (150GB! - 44,808 cells )
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63473&format=file
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63473/suppl/GSE63473\_RAW.tar

Didn't seem to be any actual RAW data here.

----------------------------------------------------------------------

##	20171228

###	Need to do / Doing / Did ...

```
create_ec2_instance.bash --profile syryu --key ~/.aws/JakeSYRyu.pem --instance-type m4.2xlarge

sudo yum update
sudo yum install htop

cd
git clone https://github.com/unreno/syryu/
cd ~/syryu
git pull
ln -s Makefile.example Makefile
make install
```

link all DropSeq scripts, (as well as jar dir and 3rdparty dir), to ~/.local/bin/

link STAR to ~/.local/bin/

include new sequence in mm10. New name mm10a
```
chmod 444 ~/mm10/*
mkdir ~/mm10a/
cd ~/mm10a
cp ~/mm10/mm10.fasta mm10a.fasta
chmod +w mm10a.fasta
cat ~/syryu/misc/EGFP_SV40_polya.fasta >> mm10a.fasta
chmod -w mm10a.fasta
ln -s ~/mm10/mm10.refFlat mm10a.refFlat
```

And create dict using Pcard’s CreateSequenceDictionary (page 3 in drop-seq cook book)? (Not sure why)

```
java -jar ~/picard.jar CreateSequenceDictionary R=mm10a.fasta
chmod -w mm10a.dict
```

update R and install libraries needed for Seurat

```
R
install.packages("devtools")
library(devtools)
install.packages("httpuv")
install.packages("igraph")
install_github("igraph/rigraph")
install.packages("Seurat")
```

Had to manually edit "sudo vi /usr/lib64/R/etc/Makeconf"
	and copy 4 lines like CXX1X\* to CXX11\*

Still didn't work.

`cd /usr/lib64/R/etc/; sudo cp Makeconf.rpmnew Makeconf`


Create New AMI with syryu script installed

`sudo halt`



Create new mm10a.gtf???? How? Needed?


```
create_ec2_instance.bash --profile syryu --key ~/.aws/JakeSYRyu.pem --instance-type m4.2xlarge --volume-size 100


ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[0].Instances[0].PublicIpAddress' --instance-ids i-07b8dd50baee3a6cd | tr -d '"' )
echo $ip
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
```


create STAR reference for mm10a (with mm10 gtf, not sure if will work)
```
cd ~/working/
mkdir -p ~/working/mm10a_star
STAR --genomeFastaFiles ~/mm10a/mm10a.fasta  --runMode genomeGenerate --genomeDir ~/working/mm10a_star --runThreadN 8 --sjdbGTFfile ~/mm10/mm10.gtf
mkdir -p ~/working/dropseq
mv ~/working/Log.out ~/working/dropseq/star_mm10a_creation.log
```


###	Upload data

`scp -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no  /Users/jakewendt/BaseSpace/Minkyung_1763-56931876/FASTQ*/*/*.bam ec2-user@$ip:working/`


###	RUN DATA

```
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
drop_seq.bash ~/working/*bam > drop_seq.log 2>&1 &
```


###	Download results

`rsync --archive --verbose --compress --rsh "ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem" --progress --delete ec2-user@$ip:working/dropseq/ ~/syryu/20171228a.drop_seq_alignment/`


Several memory failures in R.

----------------------------------------------------------------------

##	20171229a ...

```
create_ec2_instance.bash --profile syryu --key ~/.aws/JakeSYRyu.pem

aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress'

ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
```

Add links to other mm10 files or simply rename and include them (renamed them all)

Update git repo and install updated scripts

Create new ami


----------------------------------------------------------------------

##	20171229b ...

New processing machine with at least 10GB more memory (>45GB) (r4.2xlarge)

```
m4.2xlarge   8   26   32   EBS Only      $0.4 per Hour (NEED MORE MEMORY)
m5.4xlarge  16   61   64   EBS Only      $0.768 per Hour
r3.2xlarge   8   26   61   1 x 160 SSD   $0.665 per Hour
r4.2xlarge   8   27   61   EBS Only      $0.532 per Hour (TRY ME!)
x1e.xlarge   4   12  122   1 x 120 SSD   $0.834 per Hour
```

( and less processors? Change the STAR ref creation if number changes.)
(	Can't seem to get more memory AND less vCPU )
Any way to parallelize? Could run 2 at a time? Could run out of memory though if want it at the same time.


```
create_ec2_instance.bash --profile syryu --key ~/.aws/JakeSYRyu.pem --instance-type r4.2xlarge --volume-size 100


ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
```





###	create STAR reference for mm10a (with mm10 gtf and other files) (takes about 40 minutes)

```
cd ~/working/
mkdir -p ~/working/mm10a_star
STAR --genomeFastaFiles ~/mm10a/mm10a.fasta  --runMode genomeGenerate --genomeDir ~/working/mm10a_star --runThreadN 8 --sjdbGTFfile ~/mm10a/mm10a.gtf

mkdir -p ~/working/dropseq
mv ~/working/Log.out ~/working/dropseq/star_mm10a_creation.log
chmod 444 ~/working/mm10a_star/*
chmod 444 ~/working/dropseq/star_mm10a_creation.log
```


###	Upload data (takes about 5 minutes)

`scp -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no  /Users/jakewendt/BaseSpace/Minkyung_1763-56931876/FASTQ*/*/*.bam ec2-user@$ip:working/`


###	RUN DATA (takes about 8 hours)

```
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
cd ~/working/dropseq
drop_seq.bash ~/working/*bam > drop_seq.log 2>&1 &
```



----------------------------------------------------------------------

##	20171229c ... FULL STOP

###	New updates

Merge all the 1\*.bams into 1.bam AND SORT BY NAME
Merge all the 2\*.bams into 2.bam AND SORT BY NAME

```
samtools merge ~/singlecell/1.merged.bam */1*bam
samtools merge ~/singlecell/2.merged.bam */2*bam
```

Samtools sort -n will sort numbers numerically. Apparently DropSeq does not approve.
```
samtools sort -n -@ 3 -o 1.bam 1.merged.bam 
samtools sort -n -@ 3 -o 2.bam 2.merged.bam 
```

Use picard so can correctly incorrectly sort

```
cd ~/working/
java -jar ~/picard.jar MergeSamFiles \
	INPUT=1A_S4_L001.bam \
	INPUT=1A_S4_L002.bam \
	INPUT=1A_S4_L003.bam \
	INPUT=1A_S4_L004.bam \
	INPUT=1B_S3_L001.bam \
	INPUT=1B_S3_L002.bam \
	INPUT=1B_S3_L003.bam \
	INPUT=1B_S3_L004.bam \
	ASSUME_SORTED=true \
	SORT_ORDER=queryname \
	OUTPUT=1.bam

java -jar ~/picard.jar MergeSamFiles \
	INPUT=2A_S2_L001.bam \
	INPUT=2A_S2_L002.bam \
	INPUT=2A_S2_L003.bam \
	INPUT=2A_S2_L004.bam \
	INPUT=2B_S1_L001.bam \
	INPUT=2B_S1_L002.bam \
	INPUT=2B_S1_L003.bam \
	INPUT=2B_S1_L004.bam \
	ASSUME_SORTED=true \
	SORT_ORDER=queryname \
	OUTPUT=2.bam
```

Add 
	MIN\_NUM\_GENES\_PER\_CELL=100
	NUM\_CORE\_BARCODES=100
in DigitalExpression stage

Use "ds <- CreateSeuratObject(raw.data = ds.data, min.cells = 3,  min.genes = 200, is.expr=1)"
in seurat.R
 
```
cd ~/working/dropseq
drop_seq.bash ~/working/?.bam > drop_seq.log 2>&1 &
```

I am expecting memory failure. The merging of 8 data sets has resulted in a 1GB bam file.
While the dge is only 17MB, there is 59GB of 61GB memory used and it is still loading.

Failed

```
[1] "Loading data from error_detected.dge.txt.gz"
Error: cannot allocate vector of size 18.1 Gb
Execution halted
Warning message:
system call failed: Cannot allocate memory 
```



Perhaps try x1e.xlarge. Only 4 CPUs so may need to modify STAR ref creation all


###	Download results

`rsync --archive --verbose --compress --rsh "ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress --delete ec2-user@$ip:working/dropseq/ ~/syryu/20171229a.drop_seq_alignment/`


----------------------------------------------------------------------

##	20171229d Another memory failure. No new AMI. Just update existing


###	New processing machine with at least 100GB memory

```
m4.2xlarge   8   26   32   EBS Only      $0.4 per Hour (NEED MORE MEMORY)
r4.2xlarge   8   27   61   EBS Only      $0.532 per Hour (NEED MORE)
x1e.xlarge   4   12  122   1 x 120 SSD   $0.834 per Hour (TRY ME!)
```

Not sure if 120 SSD is only available? Max? Min? Seems irrelevant.


```
create_ec2_instance.bash --profile syryu --key ~/.aws/JakeSYRyu.pem --instance-type x1e.xlarge --volume-size 100

ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
```


###	create STAR reference for mm10a (with mm10 gtf and other files) (takes about 40-90 minutes depending on core count)

```
cd ~/working/
mkdir -p ~/working/mm10a_star
STAR --genomeFastaFiles ~/mm10a/mm10a.fasta  --runMode genomeGenerate --genomeDir ~/working/mm10a_star --runThreadN 4 --sjdbGTFfile ~/mm10a/mm10a.gtf


mkdir -p ~/working/dropseq
mv ~/working/Log.out ~/working/dropseq/star_mm10a_creation.log
chmod 444 ~/working/mm10a_star/*
chmod 444 ~/working/dropseq/star_mm10a_creation.log
```




###	Upload data (takes about 5 minutes)

`scp -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no  /Users/jakewendt/BaseSpace/Minkyung_1763-56931876/FASTQ*/*/*.bam ec2-user@$ip:working/`


###	Merge (This takes about an hour)

Merge all the 1\*.bams into 1.bam AND SORT BY NAME
Merge all the 2\*.bams into 2.bam AND SORT BY NAME

```
cd ~/working/
java -jar ~/picard.jar MergeSamFiles \
	INPUT=1A_S4_L001.bam \
	INPUT=1A_S4_L002.bam \
	INPUT=1A_S4_L003.bam \
	INPUT=1A_S4_L004.bam \
	INPUT=1B_S3_L001.bam \
	INPUT=1B_S3_L002.bam \
	INPUT=1B_S3_L003.bam \
	INPUT=1B_S3_L004.bam \
	ASSUME_SORTED=true \
	SORT_ORDER=queryname \
	OUTPUT=1.bam

java -jar ~/picard.jar MergeSamFiles \
	INPUT=2A_S2_L001.bam \
	INPUT=2A_S2_L002.bam \
	INPUT=2A_S2_L003.bam \
	INPUT=2A_S2_L004.bam \
	INPUT=2B_S1_L001.bam \
	INPUT=2B_S1_L002.bam \
	INPUT=2B_S1_L003.bam \
	INPUT=2B_S1_L004.bam \
	ASSUME_SORTED=true \
	SORT_ORDER=queryname \
	OUTPUT=2.bam
```

(Downloaded 1.bam and 2.bam so can just upload next time if needed.)

###	Cleanup

`\rm -f *_*_*.bam` 
 
###	Update!

```
cd ~/syryu
git pull
make install
```


###	RUN DATA (takes about 8 hours)

```
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
cd ~/working/dropseq
drop_seq.bash ~/working/?.bam > drop_seq.log 2>&1 &
```

###	Download results

`rsync --archive --verbose --compress --rsh "ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress --delete ec2-user@$ip:working/dropseq/ ~/syryu/20171229b.drop_seq_alignment/`

R crashed

----------------------------------------------------------------------

##	20180102a - gonna modify and process R

```
create_ec2_instance.bash --profile syryu --key ~/.aws/JakeSYRyu.pem --instance-type x1e.xlarge --volume-size 100

ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
```



###	UPLOAD the previously downloaded data

`rsync --archive --verbose --compress --rsh "ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress ~/syryu/20171229b.drop_seq_alignment/ ec2-user@$ip:working/dropseq/`

###	Modify and run seurat.R

```
cd ~/working/dropseq/1
seurat.R > seurat.log 2>&1 &

cd ~/working/dropseq/2
seurat.R > seurat.log 2>&1 &


rsync --archive --verbose --compress --rsh "ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress --delete ec2-user@$ip:working/dropseq/ ~/syryu/20180102a.drop_seq_alignment/

```


----------------------------------------------------------------------

##	20180103a - gonna modify and process R

```
create_ec2_instance.bash --profile syryu --key ~/.aws/JakeSYRyu.pem --instance-type x1e.xlarge --volume-size 100

ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
```



###	UPLOAD the previously downloaded data

`rsync --archive --verbose --compress --rsh "ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress ~/syryu/20171229b.drop_seq_alignment/ ec2-user@$ip:working/dropseq/`


Memory is a bit of an issue so don't run these at the same time.
It may cause a collision and then both will fail.
Put this all in a script and run that.

Make sure that you run it with nohup so it doesn't timeout and crash itself if left unattended.

```
cd ~/working/dropseq/1
chmod +w .
seurat.R > seurat.log 2>&1

cd ~/working/dropseq/2
chmod +w .
seurat.R > seurat.log 2>&1
```

Check on the creation of csv files.

`rsync --archive --verbose --compress --rsh "ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress --delete ec2-user@$ip:working/dropseq/ ~/syryu/20180103a.drop_seq_alignment/`

----------------------------------------------------------------------


##	20180119


Custom fasta with matching custom gtf, dict and refFlat files.



###	Prep mm10b.fasta locally (formerly mm10a)

####	Start from scratch

`wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63472/suppl/GSE63472_mm10_reference_metadata.tar.gz`


####	Initially, the 2 additional sequences were windows files with ^M's so fix

`vi eGFP.fasta SV40polya.fasta`

####	And mm10.fasta doesn't have an EOL EOF

```
cp mm10.fasta mm10b.fasta
echo >> mm10b.fasta
```

####	Drop-seq\_alignment.sh NEEDS A SINGLE REFERENCE FASTA SO MERGE IN NEW STUFF

```
cat eGFP.fasta >> mm10b.fasta
cat SV40polya.fasta >> mm10b.fasta
chmod 444 mm10b.fasta
```

###	Prep matching mm10b.gtf

####	Rename gtf to match

```
mv mm10new.gtf mm10b.gtf
chmod 444 mm10b.gtf
```


####	Modify gtf file (If haven't done already. I have done locally now.)

The last 2 lines of this mm10b.gtf file NEED TO BE ...

```
eGFP	AddedGenes	exon	1	576	.	+	0	gene_id "eGFP"; gene_name "eGFP"; transcript_id "eGFP"; transcript_name "eGFP";
SV40polya	AddedGenes	exon	1	240	.	+	0	gene_id "SV40polya"; gene_name "SV40polya"; transcript_id "SV40polya"; transcript_name "SV40polya";
```

Or like so where the tabs have been converted to pipes (for your viewing pleasure)

```
eGFP|AddedGenes|exon|1|576|.|+|0|gene_id "eGFP"; gene_name "eGFP"; transcript_id "eGFP"; transcript_name "eGFP";
SV40polya|AddedGenes|exon|1|240|.|+|0|gene_id "SV40polya"; gene_name "SV40polya"; transcript_id "SV40polya"; transcript_name "SV40polya";
```

Seems gene\_name and transcript\_name are expected/required by ConvertToRefFlat.

I did this AFTER I create the STAR reference. Hmm. I'll likely need to redo this.



###	Start new AWS instance

```
create_ec2_instance.bash --profile syryu --key ~/.aws/JakeSYRyu.pem --instance-type x1e.xlarge --volume-size 100

ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
```

###	UPDATE!

```
sudo yum update
cd ~/syryu
git pull
make install
```

###	Destroy mm10/mm10a stuff and prep for mm10b

```
chmod -R +w ~/mm10a
/bin/rm -rf ~/mm10a
mkdir ~/mm10b
```

###	UPLOAD BAM FILES AND NEW FILES FOR MAKING REFERENCE

```
ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip

scp -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/singlecell/?.bam ec2-user@$ip:working/

scp -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/singlecell/mm10b/mm10b.gtf ~/github/unreno/syryu/singlecell/mm10b/mm10b.fasta.gz ec2-user@$ip:mm10b/
```


###	And create dict using Pcard’s CreateSequenceDictionary

This only takes a few seconds

```
cd ~/mm10b/
java -jar ~/picard.jar CreateSequenceDictionary REFERENCE=mm10b.fasta
chmod -w mm10b.dict
```

###	Create refFlat file

Not too long here either.

```
cd cd ~/mm10b/
ConvertToRefFlat ANNOTATIONS_FILE=mm10b.gtf SEQUENCE_DICTIONARY=mm10b.dict OUTPUT=mm10b.refFlat
chmod -w mm10b.refFlat
```


###	create STAR reference for mm10b with just fasta and gtf (takes about 40-90 minutes depending on core count)


```
chmod -w ~/mm10b/mm10b*
cd ~/working/
mkdir -p ~/working/mm10b_star

gunzip ~/mm10b/mm10b.fasta.gz

STAR --genomeFastaFiles ~/mm10b/mm10b.fasta --runMode genomeGenerate --genomeDir ~/working/mm10b_star --sjdbGTFfile mm10b.gtf --sjdbOverhang 100 --runThreadN 4
```

Doesn't work ... `--genomeFastaFiles <( zcat ~/mm10b/mm10b.fasta.gz )` so must gunzip 

```
mkdir -p ~/working/dropseq
mv ~/working/Log.out ~/working/dropseq/star_mm10b_creation.log
chmod 444 ~/working/mm10b_star/*
chmod 444 ~/working/dropseq/star_mm10b_creation.log
rmdir ~/working/_STARtmp
```

###	RUN DATA (takes about 8 hours)

Drop-seq\_alignment.sh takes about an hour, then dge.bash/seurat.R takes about 2 hours or so.

```
ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
cd ~/working/dropseq
nohup drop_seq.bash ~/working/?.bam > drop_seq.log 2>&1 &
```


###	DOWNLOAD

```
ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip

rsync --archive --verbose --compress --rsh "ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress --delete ec2-user@$ip:working/dropseq/ ~/github/unreno/syryu/singlecell/20180119a.drop_seq_alignment/
```














----------------------------------------------------------------------


##	20180122

Seems that our eGFP sequence was incorrect. One last time?

Custom fasta with matching custom gtf, dict and refFlat files.


###	Prep mm10c.fasta locally (formerly mm10a)

####	Start from scratch

`wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63472/suppl/GSE63472_mm10_reference_metadata.tar.gz`


####	Initially, the 2 additional sequences were windows files with ^M's so fix

`vi eGFP.fasta SV40polya.fasta`

####	And mm10.fasta doesn't have an EOL EOF

```
cp mm10.fasta mm10c.fasta
echo >> mm10c.fasta
```

####	Drop-seq\_alignment.sh NEEDS A SINGLE REFERENCE FASTA SO MERGE IN NEW STUFF

```
cat eGFP.fasta >> mm10c.fasta
cat SV40polya.fasta >> mm10c.fasta
chmod 444 mm10c.fasta
```

###	Prep matching mm10c.gtf

####	Rename gtf to match

```
mv mm10new.gtf mm10c.gtf
chmod 444 mm10c.gtf
```


####	Modify gtf file (If haven't done already. I have done locally now.)

The last 2 lines of this mm10c.gtf file NEED TO BE ...

```
eGFP	AddedGenes	exon	1	576	.	+	0	gene_id "eGFP"; gene_name "eGFP"; transcript_id "eGFP"; transcript_name "eGFP";
SV40polya	AddedGenes	exon	1	240	.	+	0	gene_id "SV40polya"; gene_name "SV40polya"; transcript_id "SV40polya"; transcript_name "SV40polya";
```

Or like so where the tabs have been converted to pipes (for your viewing pleasure)

```
eGFP|AddedGenes|exon|1|576|.|+|0|gene_id "eGFP"; gene_name "eGFP"; transcript_id "eGFP"; transcript_name "eGFP";
SV40polya|AddedGenes|exon|1|240|.|+|0|gene_id "SV40polya"; gene_name "SV40polya"; transcript_id "SV40polya"; transcript_name "SV40polya";
```

Seems gene\_name and transcript\_name are expected/required by ConvertToRefFlat.

I did this AFTER I create the STAR reference. Hmm. I'll likely need to redo this.



###	Start new AWS instance

```
create_ec2_instance.bash --profile syryu --key ~/.aws/JakeSYRyu.pem --instance-type x1e.xlarge --volume-size 100

ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
```

###	UPDATE!

```
sudo yum update
cd ~/syryu
git pull
make install
```

###	Destroy mm10/mm10a stuff and prep for mm10c

```
chmod -R +w ~/mm10a
/bin/rm -rf ~/mm10a
mkdir ~/mm10c
```

###	UPLOAD BAM FILES AND NEW FILES FOR MAKING REFERENCE

```
ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip

scp -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/singlecell/?.bam ec2-user@$ip:working/

scp -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/singlecell/mm10c/mm10c.gtf.gz ~/github/unreno/syryu/singlecell/mm10c/mm10c.fasta.gz ec2-user@$ip:mm10c/
```


```
cd ~/mm10c/
gunzip mm10c.fasta.gz
gunzip mm10c.gtf.gz
```


###	And create dict using Pcard’s CreateSequenceDictionary

This only takes a few seconds

```
cd ~/mm10c/
java -jar ~/picard.jar CreateSequenceDictionary REFERENCE=mm10c.fasta
chmod -w mm10c.dict
```

###	Create refFlat file

Not too long here either.

```
cd cd ~/mm10c/
ConvertToRefFlat ANNOTATIONS_FILE=mm10c.gtf SEQUENCE_DICTIONARY=mm10c.dict OUTPUT=mm10c.refFlat
chmod -w mm10c.refFlat
```


###	create STAR reference for mm10c with just fasta and gtf (takes about 40-90 minutes depending on core count)


```
chmod -w ~/mm10c/mm10c*
cd ~/working/
mkdir -p ~/working/mm10c_star

gunzip ~/mm10c/mm10c.fasta.gz

nohup STAR --genomeFastaFiles ~/mm10c/mm10c.fasta --runMode genomeGenerate --genomeDir ~/working/mm10c_star --sjdbGTFfile ~/mm10c/mm10c.gtf --sjdbOverhang 100 --runThreadN 4 &
```

Doesn't work ... `--genomeFastaFiles <( zcat ~/mm10c/mm10c.fasta.gz )` so must gunzip 

```
mkdir -p ~/working/dropseq
mv ~/working/Log.out ~/working/dropseq/star_mm10c_creation.log
chmod 444 ~/working/mm10c_star/*
chmod 444 ~/working/dropseq/star_mm10c_creation.log
rmdir ~/working/_STARtmp
```

###	RUN DATA (takes about 8 hours)

Drop-seq\_alignment.sh takes about an hour, then dge.bash/seurat.R takes about 2 hours or so.

```
ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
cd ~/working/dropseq
nohup drop_seq.bash ~/working/?.bam > drop_seq.log 2>&1 &
```


###	DOWNLOAD

```
ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip

rsync --archive --verbose --compress --rsh "ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress --delete ec2-user@$ip:working/dropseq/ ~/github/unreno/syryu/singlecell/20180122a.drop_seq_alignment/
```





##	Post Run Notes

Probably should've stuck with a single name for the modified reference. Say mm10ryulab.


