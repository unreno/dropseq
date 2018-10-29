#	DropSeq Processing

##	20180228


###	Preprocess locally ...

* Download mm10c.fasta.gz, mm10c.refFlat.gz, mm10c.dict.gz, mm10c.gtf.gz from Azure
* Download picard.jar `wget https://github.com/broadinstitute/picard/releases/download/2.17.0/picard.jar`
* Download fastq data from Illumina
* Convert fastq data to bams with convert\_fastq\_files\_to\_merged\_sorted\_bams.bash


###	Start new AWS instance

```BASH
create_ec2_instance.bash --profile syryu --key ~/.aws/JakeSYRyu.pem --instance-type x1e.xlarge --volume-size 100

ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
```

###	UPDATE!

```BASH
sudo yum update
cd ~/syryu
git pull
make install
```

###	Destroy mm10/mm10a stuff and prep for mm10c

```BASH
chmod -R +w ~/mm10a
/bin/rm -rf ~/mm10a
mkdir ~/mm10c
```

###	UPLOAD BAM FILES AND NEW FILES FOR MAKING REFERENCE

```BASH
ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip

scp -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/drop_seq/B?.bam ec2-user@$ip:working/

scp -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/drop_seq/mm10c.*.gz ec2-user@$ip:mm10c/
```


```BASH
cd ~/mm10c/
gunzip mm10c.fasta.gz
gunzip mm10c.gtf.gz
gunzip mm10c.refFlat.gz
gunzip mm10c.dict.gz
```

###	create STAR reference for mm10c with just fasta and gtf (takes about 40-90 minutes depending on core count)


```BASH
chmod -w ~/mm10c/mm10c*
cd ~/working/
mkdir -p ~/working/mm10c_star

nohup STAR --genomeFastaFiles ~/mm10c/mm10c.fasta --runMode genomeGenerate --genomeDir ~/working/mm10c_star --sjdbGTFfile ~/mm10c/mm10c.gtf --sjdbOverhang 100 --runThreadN 4 &
```

```BASH
mkdir -p ~/working/dropseq
mv ~/working/Log.out ~/working/dropseq/star_mm10c_creation.log
chmod 444 ~/working/mm10c_star/*
chmod 444 ~/working/dropseq/star_mm10c_creation.log
rmdir ~/working/_STARtmp
```

###	RUN DATA (takes way more than 8 hours)

Drop-seq\_alignment.sh takes about 2 hours, then dge.bash/seurat.R takes about 7 hours or so.


read.table fails, likely due to memory shortage. Attempting to better compute memory needed.

First run sample 1's dge, 
```BASH
zcat 1-error_detected.dge.txt | wc -l
# 15564
zcat 1-error_detected.dge.txt  | head -1 | awk '{print NF}'
#	311392
expr 311392 \* 15564
#	4846505088

echo 76000000000 / 4846505088 | bc -l 
#	15.68140311833713688241
```

used at least 76GB so looks like 16 bytes each.

B3 has

```BASH
zcat error_detected.dge.txt.gz | head -1 | awk '{print NF}'
#	1030357
zcat error_detected.dge.txt.gz | wc -l
#	17895
expr 1030357 \* 17895
#	18438238515
echo 18438238515 * 16  | bc -l
295011816240


echo '1246168 * 18125 * 16 / 1000000000' | bc -l
361.38872
```

So looks like I need a machine with at least 300 GB

The final object doesn't use this much memory though.



```BASH
ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
cd ~/working/dropseq
nohup drop_seq.bash ~/working/*.bam > drop_seq.log 2>&1 &
```


###	DOWNLOAD

```BASH
ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip

rsync --archive --verbose --compress --rsh "ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress --delete ec2-user@$ip:working/dropseq/ ~/github/unreno/syryu/drop_seq/20180228a.drop_seq_alignment/
```





##	Post Run Notes

Probably should've stuck with a single name for the modified reference. Say mm10ryulab.

Given the computed memory need, will need the x1e.4xlarge to run seurat.R

```
x1e.xlarge	4	12	122	1 x 120 SSD	$0.834 per Hour
x1e.2xlarge	8	23	244	1 x 240 SSD	$1.668 per Hour
x1e.4xlarge	16	47	488	1 x 480 SSD	$3.336 per Hour
```

Looks like I'm gonna process pre-merging for these samples as well as the previous ones.


