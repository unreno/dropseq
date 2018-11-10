FROM ubuntu


#	Apparently, a new "image" is created after every "step" so cleanup after yourself.


#ENV HOME /root
#	stop complaints like 
#	debconf: (Can't locate Term/ReadLine.pm 
#	with ...
ENV DEBIAN_FRONTEND noninteractive

WORKDIR /root
 
RUN apt-get update \
	&& apt-get install -y apt-utils dialog bzip2 gcc gawk \
		zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev make \
		libssl-dev libncurses5-dev zip g++ git libtbb-dev wget \
		r-base r-base-dev libhdf5-dev \
	&& apt-get clean

#	using wget in RUN rather that ADD and then RUN does save a tiny bit of size on the final image.


RUN cd ~ \
	&& git clone https://github.com/unreno/dropseq \
	&& cd dropseq \
	&& ln -s Makefile.example Makefile \
	&& make BASE_DIR="/usr/local" install 
# \
#	&& cd ~ \
#	&& /bin/rm -rf dropseq


ENV STAR_VERSION 2.5.3a
ENV STAR_URL "https://github.com/alexdobin/STAR/archive"
RUN cd ~ \
	&& wget -O STAR-${STAR_VERSION}.tar.gz ${STAR_URL}/${STAR_VERSION}.tar.gz \
	&& tar xfvz STAR-${STAR_VERSION}.tar.gz \
	&& cp STAR-${STAR_VERSION}/bin/Linux_x86_64/STAR /usr/local/bin/ \
	&& cd ~ \
	&& /bin/rm -rf STAR-${STAR_VERSION}.tar.gz STAR-${STAR_VERSION} 
#\
#	&& which STAR


ENV DROP_SEQ_VERSION 1.13
ENV DROP_SEQ_URL "https://github.com/broadinstitute/Drop-seq/releases/download"
RUN cd ~ \
	&& wget ${DROP_SEQ_URL}/v${DROP_SEQ_VERSION}/Drop-seq_tools-${DROP_SEQ_VERSION}.zip \
	&& unzip Drop-seq_tools-${DROP_SEQ_VERSION}.zip \
	&& /bin/rm -rf Drop-seq_tools-${DROP_SEQ_VERSION}.zip


# wget -O STAR-2.6.1c.tar.gz https://github.com/alexdobin/STAR/archive/2.6.1c.tar.gz
# wget https://github.com/broadinstitute/Drop-seq/releases/download/v2.0.0/Drop-seq_tools-2.0.0.zip




#	REFERENCE ACQUISTION AND PREPARATION

RUN wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63472/suppl/GSE63472_mm10_reference_metadata.tar.gz \
	&& tar xvfz GSE63472_mm10_reference_metadata.tar.gz \
	&& /bin/rm -rf GSE63472_mm10_reference_metadata.tar.gz \
	&& chmod a+w mm10/* \
	&& sed -i 's/^>chr/>/' mm10/mm10.fasta


RUN cat ~/dropseq/misc/*dict >> mm10/mm10.dict \
	&& cat ~/dropseq/misc/*fasta >> mm10/mm10.fasta \
	&& cat ~/dropseq/misc/*gtf >> mm10/mm10.gtf \
	&& cat ~/dropseq/misc/*refFlat >> mm10/mm10.refFlat \
	&& chmod 444 mm10/*



#
#	This will require LOTS of memory. Perhaps at a memory check? About 20GB for mm10
#
#	Trying 16 GB / 3 GB SWAP - Nope
#
RUN mkdir $HOME/mm10STAR \
	&& STAR --genomeFastaFiles mm10/mm10.fasta --runMode genomeGenerate \
		--genomeDir $HOME/mm10STAR --sjdbGTFfile mm10/mm10.gtf --runThreadN 40



#	java?



#	write prepare_R_environment.Rscript or something 

RUN echo "update.packages(ask=FALSE);source('https://bioconductor.org/biocLite.R');biocLite( c('devtools','Seurat','pryr','gdata','optparse') );" | R --no-save


#	Create with ...
#
#		docker build -t dropseq .
#		docker run -ti dropseq
#
#	How many threads will it use?
#	How to pass it a param? a file?
#	How to get output out to the hosting system?
#
#	Additional connection to running container
#		docker exec -it CONTAINER_ID bash
#
#	To increase resources, need to modify Daemon GUI settings.
#
#
#
#	Run docker without sudo
#		Create docker group (may exist already)
#			sudo groupadd docker
#		Add user to docker group
#			sudo usermod -aG docker $USER
#		Logout / Login
#			exit
#			ssh herv
#
#	Monitor
#		docker stats
#
#	Delete ALL container
#		docker rm $( docker ps -qa )
#
#	Delete ALL images
#		docker rmi $( docker images -q )
#
