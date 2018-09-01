#	DropSeq Processing

##	20180518

New sample prep.

Downloaded Project MK1779 fastq.gz files from Illumina for "Batch6" and "Batch7" (B6 and B7)

Copied all Batch6\*fastq.gz into B6 directory and ran convert\_fastq\_files\_to\_merged\_sorted\_bams.bash.

Copied all Batch7\*fastq.gz into B7 directory and ran convert\_fastq\_files\_to\_merged\_sorted\_bams.bash.

Uploaded 2 bam files to Azure Storage



-r--r--r--  1 jakewendt 7934561487 May 18 16:11 B6.bam
-r--r--r--  1 jakewendt 4214557648 May 18 16:01 B7.bam

B6.bam is almost 8GB which is the size of 1, 2, 3 and 4 combined.
This may not be enough, but we shall see.


From last ubuntu image, create virtual machine.
*	name : ubuntu
* username : jake
* ssh public key : cat ~/.ssh/id_rsa.pub
* resource group : ubuntu

Size E64-16s_v3 (16 vcpus, 432 GB memory)

IP Address 40.114.26.27



###	UPDATE!


```BASH
ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@40.114.26.27

sudo apt update
sudo apt full-upgrade
sudo apt autoremove
sudo reboot

ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@40.114.26.27

cd ~/syryu
git pull
make install



mkdir ~/working/dropseq

cat ~/dest-key 

azcopy --destination ~/working/B6.bam --source https://ryulab.blob.core.windows.net/ryulab/DropSeq/bams/B6.bam --source-key $( cat ~/dest-key )
azcopy --destination ~/working/B7.bam --source https://ryulab.blob.core.windows.net/ryulab/DropSeq/bams/B7.bam --source-key $( cat ~/dest-key )

```



Remotely ...

Trying with 81GB free space. B6 8GB file took 70GB but still going.
This script really needs to more regularly clean up after itself rather than waiting until it completes.

```BASH
ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@40.114.26.27
cd ~/working/dropseq
nohup drop_seq.bash ~/working/B[67].bam > drop_seq.log 2>&1 &
```

Sadly, B6 crashed in R. Likely ran out of memory even with 432GB!


`/home/jake/.local/bin/dge.bash: line 54: 22662 Killed                  create_seurat.R`

Created monster swap space. Let's see if this makes a difference.

`df -h` and noticed that /dev/sdb1 was over 800GB and not mounted?

```BASH
sudo mount /dev/sdb1 ~/tmpdrive/

sudo fallocate -l 500G ~/tmpdrive/swap
sudo mkswap ~/tmpdrive/swap
sudo chmod 600 ~/tmpdrive/	swap 
sudo swapon ~/tmpdrive/swap 

sudo swapoff ~/tmpdrive/swap 
sudo rm ~/tmpdrive/swap 
sudo umount ~/tmpdrive/
```


B7 worked great, no need for the swap space.
Rerunning B6 by hand.

```BASH
cd ~/working/dropseq/B6	

nohup create_seurat.R > create_seurat.log 2>&1 &

nohup seurat.R --redo > seurat.log 2>&1 &
```

Reprocessing is using the swap space! R sucks at memory usage. It thinks that its only using 110GB but its using 600GB.

`create_seurat.R` seems to have worked but ended with the message ...?
```
Warning message:
system call failed: Cannot allocate memory 
```
Not sure what command failed.

`seurat.R` ran without a problem!



###	DOWNLOAD LOCALLY (if desired)

```BASH
rsync --archive --verbose --compress --rsh "ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress --delete jake@40.114.26.27:working/dropseq/ ~/github/unreno/syryu/drop_seq/20180518.drop_seq_alignment/
```


###	UPLOAD TO AZURE STORAGE

Browse to portal.azure.com, Storage Account -> ryulab -> Access Keys to find a key.

Cleanup and upload data to Azure Storage and prep to save VM image ...

Remotely ...

```BASH
azcopy --verbose --source ~/working/dropseq --destination https://ryulab.blob.core.windows.net/ryulab/DropSeq/20180518.drop_seq_alignment --recursive --dest-key $( cat ~/dest-key )

sudo swapoff ~/tmpdrive/swap 
sudo rm ~/tmpdrive/swap 
sudo umount ~/tmpdrive/

sudo waagent -deprovision
```

Using the web portal GUI, save the image




##################################################


##	Post Run Notes

Probably should've stuck with a single name for the modified reference. Say mm10ryulab.


