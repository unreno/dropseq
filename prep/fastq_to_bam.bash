#!/usr/bin/env bash



set -e	#	exit if any command fails
set -u	#	Error on usage of unset variables
set -o pipefail


export DROP_SEQ_PATH=~/.local/Drop-seq_tools-2.3.0
export PICARD_PATH=${DROP_SEQ_PATH}/3rdParty/picard/

#	
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_01_L001-ds.19308f894f6049eca2776fb6bffb3e1c/1A_S4_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_01_L002-ds.09195cc4f8824315a375f66e60c3e44e/1A_S4_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_01_L003-ds.cbfee5bb075b44e8bbb63830ead3b171/1A_S4_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_01_L004-ds.b03b7b9edb7b44d98048339e6690dacb/1A_S4_L004_R1_001.fastq.gz
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_02_L001-ds.90ce505fbe314b55915975a7f71ca06d/1B_S3_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_02_L002-ds.33cf5ad03fc3400f89d9afefd78b95c6/1B_S3_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_02_L003-ds.d21d322b7f504624a18ddfc1393c4a6c/1B_S3_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_02_L004-ds.afd4de8f6c0241f087ce812a63397625/1B_S3_L004_R1_001.fastq.gz
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_03_L001-ds.4c49cc01ef364c829c93d5416407c3ae/2A_S2_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_03_L002-ds.1beb740dd43f40fe94d01c68d70e01a2/2A_S2_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_03_L003-ds.ee0be2f1073b47d083679a1b81232597/2A_S2_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_03_L004-ds.d007b8db292040c393e0ee6b28018f96/2A_S2_L004_R1_001.fastq.gz
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_04_L001-ds.dc41d906aa6a4c89abe652e338efdc06/2B_S1_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_04_L002-ds.0e154a0ebf8443dcb4cf52f14da97317/2B_S1_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_04_L003-ds.4015d14d4e084a4db41c00448fb86c96/2B_S1_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/Minkyung_1763-56931876/FASTQ_Generation_2017-12-13_15_11_00Z-67441595/1763_04_L004-ds.c1b551c677944b12b1418210dcb9ece7/2B_S1_L004_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1771-66174108/FASTQ_Generation_2018-02-21_13_00_34Z-80978995/1171_2_L001-ds.a394657148c047b998a67e62a487c50e/B4_S2_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1771-66174108/FASTQ_Generation_2018-02-21_13_00_34Z-80978995/1171_2_L002-ds.6c8ca1c2b5dc4a8da3eb6311ae3e431c/B4_S2_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1771-66174108/FASTQ_Generation_2018-02-21_13_00_34Z-80978995/1171_2_L003-ds.af746dd8797a46688cb0fe1ade4c83fe/B4_S2_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1771-66174108/FASTQ_Generation_2018-02-21_13_00_34Z-80978995/1171_2_L004-ds.9d4ec721493649a9b88c45d0f70b2613/B4_S2_L004_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1771-66174108/FASTQ_Generation_2018-02-21_13_00_34Z-80978995/1771_1_L001-ds.e7f19e2cf325454a8d62157e11d4c579/B3_S1_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1771-66174108/FASTQ_Generation_2018-02-21_13_00_34Z-80978995/1771_1_L002-ds.5e6ed93f917f4ac8bc88716001e0addf/B3_S1_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1771-66174108/FASTQ_Generation_2018-02-21_13_00_34Z-80978995/1771_1_L003-ds.662be013f2a545f6bc54e231b6d8f001/B3_S1_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1771-66174108/FASTQ_Generation_2018-02-21_13_00_34Z-80978995/1771_1_L004-ds.e910e3240b0c4fc6be35f58bc1d96877/B3_S1_L004_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1779-77564489/FASTQ_Generation_2018-05-18_14_59_47Z-96697647/1779_1_L001-ds.c3a5062301c54ecc82f4c316b4168bba/Batch6_S2_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1779-77564489/FASTQ_Generation_2018-05-18_14_59_47Z-96697647/1779_1_L002-ds.e45dee6229e34e3ea912237f87df977d/Batch6_S2_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1779-77564489/FASTQ_Generation_2018-05-18_14_59_47Z-96697647/1779_1_L003-ds.da2037da13834ac9848615b332466735/Batch6_S2_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1779-77564489/FASTQ_Generation_2018-05-18_14_59_47Z-96697647/1779_1_L004-ds.8ca7127e972d4a94bf1d95490caa33af/Batch6_S2_L004_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1779-77564489/FASTQ_Generation_2018-05-18_14_59_47Z-96697647/1779_2_L001-ds.d47e54ba85dd45448d53d2ccfc05578f/Batch7_S1_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1779-77564489/FASTQ_Generation_2018-05-18_14_59_47Z-96697647/1779_2_L002-ds.75ae3ac16f57475ba36c3623f8baedf5/Batch7_S1_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1779-77564489/FASTQ_Generation_2018-05-18_14_59_47Z-96697647/1779_2_L003-ds.4e405c3a4a914a49bcd419c415bbed93/Batch7_S1_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK_1779-77564489/FASTQ_Generation_2018-05-18_14_59_47Z-96697647/1779_2_L004-ds.e0c3a28755194459a174054e2feb3b99/Batch7_S1_L004_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-1_L001-ds.d4ce8144325748edac6780dbd8513b69/Batch8-1_S2_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-1_L002-ds.0b126abe07334522a28b12cf0fecd70c/Batch8-1_S2_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-1_L003-ds.5481bac9d78c4dbba65f4a52391e2315/Batch8-1_S2_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-1_L004-ds.7b6bbbfe40004c1e83194d27637be54c/Batch8-1_S2_L004_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-2_L001-ds.5cff5b807252422e9307b363602cd4f0/Batch8-2_S3_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-2_L002-ds.82cb7486f0174906b55c96e3f183678e/Batch8-2_S3_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-2_L003-ds.0b52f2937eef44f4b7760e572baa2bc0/Batch8-2_S3_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-2_L004-ds.6e1acb4525f349aab9a79774fa9b0805/Batch8-2_S3_L004_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-3_L001-ds.749fed820c3c45adbdd866e035d9be20/Batch8-3_S4_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-3_L002-ds.2f694c32bdec4397af573be6da0b0057/Batch8-3_S4_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-3_L003-ds.8baa9275c9574172961c36e6af50e13b/Batch8-3_S4_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-3_L004-ds.ae57a219dfe742b39073ad0eca60b5c1/Batch8-3_S4_L004_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-4_L001-ds.bf925f11db87498f9b76c2fc5d608fcb/Batch8-4_S6_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-4_L002-ds.42feb988a07142829b9124ab9a16b9d3/Batch8-4_S6_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-4_L003-ds.31141adaff194a089a2167b5bfccad42/Batch8-4_S6_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-4_L004-ds.9eee6a40ec0643e9ab659e8d662356dd/Batch8-4_S6_L004_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-5_L001-ds.3bd8d6477923457fbd6f797f0ab20646/Batch9-1_S8_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-5_L002-ds.66796291cd15410492f3842d9cdd0131/Batch9-1_S8_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-5_L003-ds.595f8891639b48d09b2ef34e0c0ec346/Batch9-1_S8_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-5_L004-ds.8a81ff2cd4ae499c8fb74c221a0a5632/Batch9-1_S8_L004_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-6_L001-ds.3f9048360dc44e6bb3a8cf130952667c/Batch9-2_S1_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-6_L002-ds.cf91ac8daf9d48f1abfc5637f83a787f/Batch9-2_S1_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-6_L003-ds.dd1cc15bb7e04cf7ad952cd5092cc095/Batch9-2_S1_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-6_L004-ds.4a34bb563ced42ffb2cafb2a4ab6f3c2/Batch9-2_S1_L004_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-7_L001-ds.82756f529a4043219403fde2907af13b/Batch9-3_S7_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-7_L002-ds.498b0271ff244c2491ae657fd2768464/Batch9-3_S7_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-7_L003-ds.87cbdd4866894256a8f789c1cfd845d5/Batch9-3_S7_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-7_L004-ds.2413c589dce544f59ce7de25be0595d5/Batch9-3_S7_L004_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-8_L001-ds.b96d081df65f4d698ae545aafb05aed7/Batch9-4_S5_L001_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-8_L002-ds.9b736e1837a449318797527d83cd3744/Batch9-4_S5_L002_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-8_L003-ds.e2a1310643244a47840483e2ac7fcb4a/Batch9-4_S5_L003_R1_001.fastq.gz
#	/home/jake/BaseSpace/MK1803-108592487/FASTQ_Generation_2018-12-13_13_11_26Z-143234098/1803-8_L004-ds.24cc66b155ff4164a908af78ac52b49a/Batch9-4_S5_L004_R1_001.fastq.gz
#	

for r1 in ~/BaseSpace/*/FASTQ*/*/*_R1_001.fastq.gz ; do
	r2=${r1/_R1_/_R2_}

	base=$( basename $r1 _R1_001.fastq.gz )
	base=B${base}
	base=${base/BB/B}
	base=${base/atch/}

	f=${base}.bam
	if [ -f $f ] && [ ! -w $f ] ; then
		echo "Write-protected $f exists. Skipping."
	else
		echo "Creating $f"
		java -jar $PICARD_PATH/picard.jar FastqToSam \
			F1=${r1} \
			F2=${r2} \
			O=${f} SM=${base}
		chmod a-w $f
	fi

done

for l1 in *_L001.bam ; do
	l1=${l1/_S?/}
	l2=${l1/_L001/_L002}
	l3=${l1/_L001/_L003}
	l4=${l1/_L001/_L004}

	base=$( basename $l1 _L001.bam )

	f=${base}.bam
	if [ -f $f ] && [ ! -w $f ] ; then
		echo "Write-protected $f exists. Skipping."
	else
		echo "Creating $f"
		java -jar $PICARD_PATH/picard.jar MergeSamFiles \
			INPUT=${l1} \
			INPUT=${l2} \
			INPUT=${l3} \
			INPUT=${l4} \
			ASSUME_SORTED=true \
			SORT_ORDER=queryname \
			OUTPUT=${base}.bam
		chmod a-w $f
	fi

done


for a in B?A.bam ; do
	echo $a
	base=$( basename $a A.bam )

	f=${base}.bam
	if [ -f $f ] && [ ! -w $f ] ; then
		echo "Write-protected $f exists. Skipping."
	else
		echo "Creating $f"
		java -jar $PICARD_PATH/picard.jar MergeSamFiles \
			INPUT=${base}A.bam \
			INPUT=${base}B.bam \
			ASSUME_SORTED=true \
			SORT_ORDER=queryname \
			OUTPUT=${base}.bam
		chmod a-w $f
	fi

done

for a in B?-1.bam ; do
	echo $a
	base=$( basename $a -1.bam )

	f=${base}.bam
	if [ -f $f ] && [ ! -w $f ] ; then
		echo "Write-protected $f exists. Skipping."
	else
		echo "Creating $f"
		java -jar $PICARD_PATH/picard.jar MergeSamFiles \
			INPUT=${base}-1.bam \
			INPUT=${base}-2.bam \
			INPUT=${base}-3.bam \
			INPUT=${base}-4.bam \
			ASSUME_SORTED=true \
			SORT_ORDER=queryname \
			OUTPUT=${base}.bam
		chmod a-w $f
	fi

done


for bam in B?.bam ; do
	echo $bam
	base=$( basename $bam .bam )

	f=${base}
	if [ -d $f ] && [ ! -w $f ] ; then
		echo "Write-protected $f exists. Skipping."		#	DIRECTORY
	else
		echo "Creating $f"
#		drop_seq.bash --drop_seq ${DROP_SEQ_PATH}/Drop-seq_alignment.sh \
#			--star ~/.local/STAR-2.7.1a/bin/Linux_x86_64/STAR \
#			--genomedir ~/.github/unreno/dropseq/prep/mm10x.STAR-2.7.1a \
#			--referencefasta ~/.github/unreno/dropseq/prep/mm10x/mm10x.fasta \
#			${base}.bam > ${base}.log
#			#--estimated-num-cells 20000 \
#		chmod -R a-w $f
	fi

done
