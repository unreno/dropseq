#	DropSeq Processing

##	20180302


###	Preprocess locally ...

* Download mm10c.fasta.gz, mm10c.refFlat.gz, mm10c.dict.gz, mm10c.gtf.gz from Azure
* Download picard.jar `wget https://github.com/broadinstitute/picard/releases/download/2.17.0/picard.jar`
* Download fastq data from Illumina
* Convert B3 and B4 fastq data to bams
* Merge some of the previous sample data

```BASH
java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=1_L001.bam INPUT=1A_S4_L001.bam INPUT=1B_S3_L001.bam

java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=1_L002.bam INPUT=1A_S4_L002.bam INPUT=1B_S3_L002.bam

java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=1_L003.bam INPUT=1A_S4_L003.bam INPUT=1B_S3_L003.bam

java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=1_L004.bam INPUT=1A_S4_L004.bam INPUT=1B_S3_L004.bam

java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=2_L001.bam INPUT=2A_S2_L001.bam INPUT=2B_S1_L001.bam

java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=2_L002.bam INPUT=2A_S2_L002.bam INPUT=2B_S1_L002.bam

java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=2_L003.bam INPUT=2A_S2_L003.bam INPUT=2B_S1_L003.bam

java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=2_L004.bam INPUT=2A_S2_L004.bam INPUT=2B_S1_L004.bam
```



###	Gonna have a go at a Ubuntu Server on Azure


Started Standard_E16s_v3, 16CPU 128GB memory, 128GB disk via the web site.




###	UPDATE!

```BASH
ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@52.168.35.1
sudo apt update
sudo apt full-upgrade
sudo apt install make

mkdir -p ~/working
mkdir -p ~/mm10c


git clone https://github.com/unreno/syryu
cd ~/syryu
ln -s Makefile.example Makefile
make install
```

```BASH
scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/drop_seq/bams/*.bam jake@52.168.35.1:working/

scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/drop_seq/mm10c.*.gz jake@52.168.35.1:mm10c/
```


```BASH

ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@52.168.35.1

sudo apt install openjdk-8-jdk

sudo apt install unzip




wget -O Drop-seq_tools.zip http://mccarrolllab.com/download/1276/
unzip Drop-seq_tools.zip
/bin/rm -rf Drop-seq_tools.zip

wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz
tar xvfz 2.5.3a.tar.gz
/bin/rm -rf 2.5.3a.tar.gz 

wget https://github.com/broadinstitute/picard/releases/download/2.17.0/picard.jar

mkdir ~/bin
cp STAR-2.5.3a/bin/Linux_x86_64_static/STAR bin/

sudo apt install htop
sudo apt install openssl libssl-dev libcurl4-openssl-dev libssh2-1-dev


mkdir .R
cat > .Renviron <<EOF
R_LIBS="/home/jake/.R"
R_LIBS_USER="/home/jake/.R"
EOF

cat > .Rprofile <<EOF
local({
	r <- getOption("repos")
  r["CRAN"] <- "https://cloud.r-project.org/"
	options(repos = r)
})
EOF


#sudo apt install r-base-core
#	apt version of R is 3.2.?. Too old
wget https://cloud.r-project.org/src/base/R-3/R-3.4.3.tar.gz

sudo apt install libx11-dev xorg-dev
cd R-3.4.3
./configure --prefix ~/.local
make install



R
install.packages("devtools")
library(devtools)
install.packages("httpuv")
#install.packages("igraph")
install_github("igraph/rigraph")
install.packages("Seurat")
install.packages("pryr")



cd ~/mm10c/
gunzip mm10c.fasta.gz
gunzip mm10c.gtf.gz
gunzip mm10c.refFlat.gz
gunzip mm10c.dict.gz


chmod -w ~/mm10c/mm10c*
cd ~/working/
mkdir -p ~/working/mm10c_star

nohup STAR --genomeFastaFiles ~/mm10c/mm10c.fasta --runMode genomeGenerate --genomeDir ~/working/mm10c_star --sjdbGTFfile ~/mm10c/mm10c.gtf --sjdbOverhang 100 --runThreadN 16 &

```








```BASH
mkdir -p ~/working/dropseq
mv ~/working/Log.out ~/working/dropseq/star_mm10c_creation.log
chmod 444 ~/working/mm10c_star/*
chmod 444 ~/working/dropseq/star_mm10c_creation.log
rmdir ~/working/_STARtmp
```

###	RUN DATA (individually, should take under an hour on Azure 16/128)

Drop-seq\_alignment.sh takes about 20 min, then dge.bash/seurat.R takes about 20 min or so.


```BASH
ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@52.168.35.1
cd ~/working/dropseq
nohup drop_seq.bash ~/working/*.bam > drop_seq.log 2>&1 &
```


###	DOWNLOAD

```BASH
rsync --archive --verbose --compress --rsh "ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress --delete jake@52.168.35.1:working/dropseq/ ~/github/unreno/syryu/drop_seq/20180302a.drop_seq_alignment/
```


##################################################


##	Post Run Notes

Probably should've stuck with a single name for the modified reference. Say mm10ryulab.


```
rename B4_S2 B4 B4_S2_L00*
rename B3_S1 B3 B3_S1_L00*
```


Sadly, failed on 6 of the 16.

Restarting on an 8 CPU 256GB

Standard_E32-8s_v3 (8 vcpus, 256 GB memory)

New ip address 13.92.226.126



```
B4\* samples ran out of memory.
system call failed: Cannot allocate memory 

The 2 B3 samples failed differently. Not sure if these are memory issues or not.
[1] "PCHeatmap"
Error in sx[1:num.genes, , drop = FALSE] : subscript out of bounds
Calls: PCHeatmap ... DimTopCells -> unique -> unlist -> lapply -> FUN -> rownames
Execution halted

[1] "PCHeatmap"
Error in sx[1:num.genes, , drop = FALSE] : subscript out of bounds
Calls: PCHeatmap ... DimTopCells -> unique -> unlist -> lapply -> FUN -> rownames
Execution halted
Warning message:
system call failed: Cannot allocate memory 
```




```BASH
ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@13.92.226.126
cd ~/working/dropseq
nohup drop_seq.bash ~/working/{B4_L00?,B3_L001,B3_L004}.bam > drop_seq.2.log 2>&1 &
```

If all 6 complete successfully, ...

All 6 create seurat object successfully, however 5 of the 6 failed at other points.


```BASH
cd ~/working/dropseq
nohup seurat_group.bash 1 2 B3 B4 > seurat_group.log 2>&1 &
```








###	DOWNLOAD

```BASH
rsync --archive --verbose --compress --rsh "ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress --delete jake@13.92.226.126:working/dropseq/ ~/github/unreno/syryu/drop_seq/20180302a.drop_seq_alignment/
```


### One more time

Gonna start a monster and see if can process the whole b3 and b4 samples

Standard_E64-16s_v3 (16 vcpus, 432 GB memory)


52.226.129.154

mkdir -p ~/working/drop_seq/B3
mkdir -p ~/working/drop_seq/B4

```BASH
rsync --archive --verbose --compress --rsh "ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress ~/github/unreno/syryu/drop_seq/20180228a.drop_seq_alignment jake@52.226.129.154:working/

scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/drop_seq/20180228a.drop_seq_alignment/B3/error_detected.dge.txt.gz jake@52.226.129.154:working/drop_seq/B3
scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/drop_seq/20180228a.drop_seq_alignment/B4/error_detected.dge.txt.gz jake@52.226.129.154:working/drop_seq/B4
```


For B4, initial seurat.R read.table and creation of seurat object took about 8 hours and used 384GB of memory


```BASH
ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@52.226.129.154
cd ~/working/20180228a.drop_seq_alignment/B3
nohup seurat.R > seurat.log 2>&1 &


cd ~/working/20180228a.drop_seq_alignment/B4
nohup seurat.R > seurat.log 2>&1 &
nohup seurat.R --redo > seurat.log 2>&1 &
```


###	DOWNLOAD

```BASH
rsync --archive --verbose --compress --rsh "ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress --delete jake@52.226.129.154:working/dropseq/ ~/github/unreno/syryu/drop_seq/20180304a.drop_seq_alignment/
```





### Post run notes

It seems that R doesn't clean up well.
Even after deleting object, it still takes up quite a bit of memory.
While the memory used doesn't show in mem\_used(), the system shows it as mostly taken.
Later in the script, another function can run out of memory and the script crashes out.
Rerunning seurat.R using my option --redo, loads the data from the stored data and never loads the dge file theu never fills the memory.
Success.
Perhaps future runs should have one script simply convert the dge to seurat and then quit leaving the analysis to load the saved seurat object and analyze.


