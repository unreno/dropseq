## CYGWIN Notes

Apparently, this was not designed to run on Cygwin.

I have had to modify several system files and the STAR source to get it to compile.

The Drop Seq script all required `cygpath` modifications.

Several of the R packages just don't.




###	Compiling STAR on cygwin

STAR 2.5.3a has many complications ...


```
SharedMemory.cpp:109:59: error: ‘SHM_NORESERVE’ was not declared in this scope

SharedMemory.cpp:254:72: error: ‘SHM_NORESERVE’ was not declared in this scope

SharedMemory.cpp

Commenting out these if/endif lines

//#ifdef COMPILE_FOR_MAC
  #define SHM_NORESERVE 0
//#endif
```







```
Parameters_openReadsFiles.cpp:82:23: error: ‘vfork’ was not declared in this scope

Change vfork to fork

Parameters_openReadsFiles.cpp

//pid_t PID=vfork();
pid_t PID=fork();
```




```
Parameters_closeReadsFiles.cpp:10:3: error: ‘kill’ was not declared in this scope

Comment out the if/endif around kill when compiling STAR in /usr/include/sys/signal.h

//#if __POSIX_VISIBLE
int kill (pid_t, int);
//#endif
```


This is moot now...


At some point, I get "EXITING because of FATAL ERROR in reads input: quality string length is not equal to sequence length" when it just isn't true so commenting out the whole block.


/*

        if ((uint) readInStream.gcount() != LreadOriginal+1) {//inconsistent read sequence and quality
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in reads input: quality string length is not equal to sequence length\n";
            errOut << readName<<"\n";
            errOut << Seq <<"\n";
            errOut << Qual <<"\n";
            errOut << "SOLUTION: fix your fastq file\n";
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
        };
*/


Seems that the input fastq quality string can commonly be terminated with a ControlM which is an additional character. Hmmm.



Adding the sed in the Drop Seq script removes the need for the above modification and concern.





###	Modify Drop Seq scripts


**Be advised:**

There a a number of modifications that need made to run this on cygwin linux on a Windows machine.

Basically, cygwin is a Linux environment, but java requires Windows paths.


```
All java wrappers
< jar_deploy_dir=$thisdir/jar
---
> jar_deploy_dir=$( cygpath -w $thisdir/jar )
```



```
Drop-seq_alignment.sh
> #picard_jar=${dropseq_root}/3rdParty/picard/picard.jar
> picard_jar=$(cygpath -w ${dropseq_root}/3rdParty/picard/picard.jar )


Add ... reference=$( cygpath -w $reference )
before its use



FASTQ Quality lines end up with trailing ^M??? (This is Control-V Carriage-Return)

sed -i -e 's/^M//' $tmpdir/unaligned_mc_tagged_polyA_filtered.fastq


```




###	Preparing R and it's packages


Cygwin includes an R library, but can't seem to install any packages.

This requires the additional cygwin developer packages of x11, pcre, readline, bz2, lzma, curl, etc. for configure/make.




  CC       src/unix/libuv_la-bsd-ifaddrs.lo
src/unix/bsd-ifaddrs.c: In function 'uv__ifaddr_exclude':
src/unix/bsd-ifaddrs.c:45:41: error: 'AF_LINK' undeclared (first use in this function); did you mean 'AF_HYLINK'?
     return (ent->ifa_addr->sa_family != AF_LINK);
                                         ^~~~~~~
                                         AF_HYLINK
src/unix/bsd-ifaddrs.c:45:41: note: each undeclared identifier is reported only once for each function it appears in
make[1]: *** [Makefile:2258: src/unix/libuv_la-bsd-ifaddrs.lo] Error 1
make[1]: Leaving directory '/tmp/RtmpHGTOPU/R.INSTALL251855f5d24f/fs/src/libuv'
make: *** [Makevars:38: libuv/.libs/libuv.a] Error 2
ERROR: compilation failed for package ‘fs’
* removing ‘/usr/lib/R/site-library/fs’
* installing *source* package ‘plyr’ ...
** package ‘plyr’ successfully unpacked and MD5 sums checked


Seems the heart of the problem is httpuv ...

Untarred httpuv 
Manual mod of httpuv
https://github.com/libuv/libuv/commit/63de1ecad3252d3e9ed2fe960c21d9387615fa45
Basically delete around line 45 in src/unix/bsd-ifaddrs.c ...
 if (exclude_type == UV__EXCLUDE_IFPHYS)
    return (ent->ifa_addr->sa_family != AF_LINK);	


re-tarred
then 
install.packages("/tmp/Rtmp5OS6GQ/downloaded_packages/httpuv_1.4.5.tar.gz")








The "fs" package contains its own copy of this code and its problems


  CC       src/unix/libuv_la-bsd-ifaddrs.lo
src/unix/bsd-ifaddrs.c: In function 'uv__ifaddr_exclude':
src/unix/bsd-ifaddrs.c:45:41: error: 'AF_LINK' undeclared (first use in this function); did you mean 'AF_HYLINK'?
     return (ent->ifa_addr->sa_family != AF_LINK);
                                         ^~~~~~~
                                         AF_HYLINK
src/unix/bsd-ifaddrs.c:45:41: note: each undeclared identifier is reported only once for each function it appears in
make[1]: *** [Makefile:2258: src/unix/libuv_la-bsd-ifaddrs.lo] Error 1
make[1]: Leaving directory '/tmp/RtmpcxemJp/R.INSTALL1f047e45da5/fs/src/libuv'
make: *** [Makevars:38: libuv/.libs/libuv.a] Error 2
ERROR: compilation failed for package ‘fs’
* removing ‘/usr/lib/R/site-library/fs’


fs/src/libuv/src/unix/bsd-ifaddrs.c

Deleting lines ..
  if (exclude_type == UV__EXCLUDE_IFPHYS)
    return (ent->ifa_addr->sa_family != AF_LINK);

> install.packages("/tmp/Rtmp5OS6GQ/downloaded_packages/fs_1.2.6.tar.gz")





And now 




  CCLD     libuv.la
make[1]: Leaving directory '/tmp/RtmpIMHVQe/R.INSTALL2b606a931bea/fs/src/libuv'
g++ -shared -L/usr/lib/R/lib -o fs.dll error.o dir.o fs.o link.o path.o file.o utils.o id.o unix/getmode.o RcppExports.o ./libuv/.libs/libuv.a -L/usr/lib/R/lib -lR -lintl -lpcre -llzma -lbz2 -lz -ltirpc -lrt -ldl -lm -liconv -licuuc -licui18n
unix/getmode.o: In function `getmode_(char const*, unsigned int)':
/tmp/RtmpIMHVQe/R.INSTALL2b606a931bea/fs/src/unix/getmode.cc:16: undefined reference to `setmode'
/tmp/RtmpIMHVQe/R.INSTALL2b606a931bea/fs/src/unix/getmode.cc:16:(.text+0xd): relocation truncated to fit: R_X86_64_PC32 against undefined symbol `setmode'
unix/getmode.o: In function `strmode_(unsigned int)':
/tmp/RtmpIMHVQe/R.INSTALL2b606a931bea/fs/src/unix/getmode.cc:28: undefined reference to `strmode'
/tmp/RtmpIMHVQe/R.INSTALL2b606a931bea/fs/src/unix/getmode.cc:28:(.text+0x76): relocation truncated to fit: R_X86_64_PC32 against undefined symbol `strmode'
collect2: error: ld returned 1 exit status
make: *** [/usr/lib/R/share/make/shlib.mk:6: fs.dll] Error 1
ERROR: compilation failed for package ‘fs’
* removing ‘/usr/lib/R/site-library/fs’
Warning message:
In install.packages("/tmp/Rtmp5OS6GQ/downloaded_packages/fs_1.2.6.tar.gz") :
  installation of package ‘/tmp/Rtmp5OS6GQ/downloaded_packages/fs_1.2.6.tar.gz’ had non-zero exit status
>



fs/src/unix/getmode.cc

//#if (defined(__APPLE__) && defined(__MACH__)) || defined(__BSD__)
//#include <string.h> /* for strmode */
//#include <unistd.h> /* for getmode / setmode */
//#else
#include "bsd/string.h" /* for strmode */
#include "bsd/unistd.h" /* for getmode / setmode */
//#endif

Try again

> install.packages("/tmp/Rtmp5OS6GQ/downloaded_packages/fs_1.2.6.tar.gz")








###	Creating data sets need to use proper paths as well

```BASH
export PICARD_PATH=~/Drop-seq_tools-1.13/3rdparty/picard/
java -jar picard.jar FastqToSam \
    F1=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_01-78286209/1A_S4_L001_R1_001.fastq.gz ) \
    F2=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_01-78286209/1A_S4_L001_R2_001.fastq.gz ) \
    SM=B1A_S4_L001 O=B1A_S4_L001.bam
java -jar picard.jar FastqToSam \
    F1=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_01-78286209/1A_S4_L002_R1_001.fastq.gz ) \
    F2=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_01-78286209/1A_S4_L002_R2_001.fastq.gz ) \
    SM=B1A_S4_L002 O=B1A_S4_L002.bam
java -jar picard.jar FastqToSam \
    F1=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_01-78286209/1A_S4_L003_R1_001.fastq.gz ) \
    F2=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_01-78286209/1A_S4_L003_R2_001.fastq.gz ) \
    SM=B1A_S4_L003 O=B1A_S4_L003.bam
java -jar picard.jar FastqToSam \
    F1=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_01-78286209/1A_S4_L004_R1_001.fastq.gz ) \
    F2=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_01-78286209/1A_S4_L004_R2_001.fastq.gz ) \
    SM=B1A_S4_L004 O=B1A_S4_L004.bam
java -jar picard.jar MergeSamFiles \
    INPUT=B1A_S4_L001.bam \
    INPUT=B1A_S4_L002.bam \
    INPUT=B1A_S4_L003.bam \
    INPUT=B1A_S4_L004.bam \
    ASSUME_SORTED=true \
    SORT_ORDER=queryname \
    OUTPUT=B1A.bam



export PICARD_PATH=~/Drop-seq_tools-1.13/3rdparty/picard/
java -jar $( cygpath -w $PICARD_PATH/picard.jar ) FastqToSam \
    F1=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_02-78283210/1B_S3_L001_R1_001.fastq.gz ) \
    F2=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_02-78283210/1B_S3_L001_R2_001.fastq.gz ) \
    SM=B1B_S3_L001 O=B1B_S3_L001.bam
java -jar $( cygpath -w $PICARD_PATH/picard.jar ) FastqToSam \
    F1=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_02-78283210/1B_S3_L002_R1_001.fastq.gz ) \
    F2=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_02-78283210/1B_S3_L002_R2_001.fastq.gz ) \
    SM=B1B_S3_L002 O=B1B_S3_L002.bam
java -jar $( cygpath -w $PICARD_PATH/picard.jar ) FastqToSam \
    F1=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_02-78283210/1B_S3_L003_R1_001.fastq.gz ) \
    F2=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_02-78283210/1B_S3_L003_R2_001.fastq.gz ) \
    SM=B1B_S3_L003 O=B1B_S3_L003.bam
java -jar $( cygpath -w $PICARD_PATH/picard.jar ) FastqToSam \
    F1=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_02-78283210/1B_S3_L004_R1_001.fastq.gz ) \
    F2=$(cygpath -w /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_02-78283210/1B_S3_L004_R2_001.fastq.gz ) \
    SM=B1B_S3_L004 O=B1B_S3_L004.bam
java -jar $( cygpath -w $PICARD_PATH/picard.jar ) MergeSamFiles \
    INPUT=B1B_S3_L001.bam \
    INPUT=B1B_S3_L002.bam \
    INPUT=B1B_S3_L003.bam \
    INPUT=B1B_S3_L004.bam \
    ASSUME_SORTED=true \
    SORT_ORDER=queryname \
    OUTPUT=B1B.bam



export PICARD_PATH=~/Drop-seq_tools-1.13/3rdparty/picard/

cd /cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_02-78286210/
java -jar $( cygpath -w $PICARD_PATH/picard.jar ) FastqToSam \
    F1=1B_S3_L001_R1_001.fastq.gz \
    F2=1B_S3_L001_R2_001.fastq.gz \
    SM=B1B_S3_L001 O=B1B_S3_L001.bam
java -jar $( cygpath -w $PICARD_PATH/picard.jar ) FastqToSam \
    F1=1B_S3_L002_R1_001.fastq.gz \
    F2=1B_S3_L002_R2_001.fastq.gz \
    SM=B1B_S3_L002 O=B1B_S3_L002.bam
java -jar $( cygpath -w $PICARD_PATH/picard.jar ) FastqToSam \
    F1=1B_S3_L003_R1_001.fastq.gz \
    F2=1B_S3_L003_R2_001.fastq.gz \
    SM=B1B_S3_L003 O=B1B_S3_L003.bam
java -jar $( cygpath -w $PICARD_PATH/picard.jar ) FastqToSam \
    F1=1B_S3_L004_R1_001.fastq.gz \
    F2=1B_S3_L004_R2_001.fastq.gz \
    SM=B1B_S3_L004 O=B1B_S3_L004.bam
java -jar $( cygpath -w $PICARD_PATH/picard.jar ) MergeSamFiles \
    INPUT=B1B_S3_L001.bam \
    INPUT=B1B_S3_L002.bam \
    INPUT=B1B_S3_L003.bam \
    INPUT=B1B_S3_L004.bam \
    ASSUME_SORTED=true \
    SORT_ORDER=queryname \
    OUTPUT=B1B.bam





export DROP_SEQ_PATH=~/Drop-seq_tools-1.13cygwin
nohup drop_seq.bash --drop_seq ${DROP_SEQ_PATH}/Drop-seq_alignment.sh \
    --estimated-num-cells 20000 --genomedir ${PWD}/mm10STAR2.5.3a \
    --referencefasta ${PWD}/mm10/mm10.fasta B1A.bam > B1A.log &





/cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_02-78283210:
total 1376640
-rwxrwx---+ 1 soyoungr None 194189543 Dec 19  2017 1B_S3_L001_R1_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 208022734 Dec 19  2017 1B_S3_L001_R2_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 153419289 Dec 19  2017 1B_S3_L002_R1_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 158434582 Dec 19  2017 1B_S3_L002_R2_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 144747122 Dec 19  2017 1B_S3_L003_R1_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 153774418 Dec 19  2017 1B_S3_L003_R2_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 191244578 Dec 19  2017 1B_S3_L004_R1_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 205539669 Dec 19  2017 1B_S3_L004_R2_001.fastq.gz

/cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_03-78269249:
total 207936
-rwxrwx---+ 1 soyoungr None 30056456 Dec 19  2017 2A_S2_L001_R1_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 32072428 Dec 19  2017 2A_S2_L001_R2_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 23106427 Dec 19  2017 2A_S2_L002_R1_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 23819248 Dec 19  2017 2A_S2_L002_R2_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 20481340 Dec 19  2017 2A_S2_L003_R1_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 21747864 Dec 19  2017 2A_S2_L003_R2_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 29564020 Dec 19  2017 2A_S2_L004_R1_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 31792433 Dec 19  2017 2A_S2_L004_R2_001.fastq.gz

/cygdrive/c/BaseSpace/Minkyung_1763-56931876/1763_04-78276228:
total 1540352
-rwxrwx---+ 1 soyoungr None 200090827 Dec 19  2017 2B_S1_L001_R1_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 213420214 Dec 19  2017 2B_S1_L001_R2_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 184441543 Dec 19  2017 2B_S1_L002_R1_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 189403446 Dec 19  2017 2B_S1_L002_R2_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 170325101 Dec 19  2017 2B_S1_L003_R1_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 179538150 Dec 19  2017 2B_S1_L003_R2_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 212254209 Dec 19  2017 2B_S1_L004_R1_001.fastq.gz
-rwxrwx---+ 1 soyoungr None 227558025 Dec 19  2017 2B_S1_L004_R2_001.fastq.gz
```








```BASH
nohup drop_seq.bash --drop_seq ${DROP_SEQ_PATH}/Drop-seq_alignment.sh --estimated-num-cells 20000 --genomedir  ${PWD}/mm10cSTAR2.5.3a --referencefasta ${PWD}/mm10c/mm10c.fasta B1B.bam > B1B.log &
```
