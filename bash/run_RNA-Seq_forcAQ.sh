#!/bin/bash

GENPATH="~/QNAP2/GENCODEv29"

echo "#################################"
echo "Prepare STAR and RSEM:"
date
echo "#################################"

STAR --runMode genomeGenerate \
        --runThreadN 18 \
        --genomeDir ${GENPATH}/STARIndex \
        --genomeFastaFiles ${GENPATH}/GRCh38.primary_assembly.genome.fa \
        --sjdbGTFfile ${GENPATH}/gencode.v29.annotation.gtf \
        --sjdbOverhang 100 


rsem-prepare-reference --gtf ${GENPATH}/gencode.v29.annotation.gtf \
       ${GENPATH}/GRCh38.primary_assembly.genome.fa \
       ${GENPATH}/RSEMref/RSEMref


ulimit -n 65536 # for sorting bam by STAR
LC_COLLATE=C # for sorting bedGraph

for LID in `cat PR2056_ID.list`
do

echo "#################################"
echo "RNA-Seq Data processing of $LID was starting at:"
date
echo "#################################"



mkdir ./OUT
echo "#################################"
echo "QualityCheking and mapping was started of $LID:"
date
echo "#################################"

rabbit_qc -w 18 -i ./RNA-Seq/RAW/${LID}_L004_R1.fastq.gz \
	-I ./RNA-Seq/RAW/${LID}_L004_R2.fastq.gz\
       	-o ./OUT/${LID}_L004_R1_tr.fastq.gz\
       	-O ./OUT/${LID}_L004_R2_tr.fastq.gz\
	-j ./OUT/${LID}_RabbitQC.json -h ./OUT/${LID}_RabbitQC.html



echo "#################################"
echo "Running mapping job by STAR:"
date
echo "#################################"

#Options imitate ENCODE's ones.
STARparCommon=" --genomeDir ${GENPATH}/STARIndex \
      	--readFilesIn ./OUT/${LID}_L004_R1_tr.fastq.gz ./OUT/${LID}_L004_R2_tr.fastq.gz \
	--outFileNamePrefix ./OUT/
	--outSAMunmapped Within --outFilterType BySJout \
	--outSAMattributes NH HI AS NM MD \
     	--outFilterMultimapNmax 20   --outFilterMismatchNmax 999 \
	--outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20 \
      	--alignIntronMax 1000000   --alignMatesGapMax 1000000   \
	--alignSJoverhangMin 8   --alignSJDBoverhangMin 1 \
	--sjdbScore 1 --readFilesCommand zcat"

# STAR parameters: run-time, controlled by DCC
STARparRun=" --runThreadN 18 --limitBAMsortRAM 60000000000"

STARparBAM="--outSAMtype BAM Unsorted --quantMode TranscriptomeSAM"
STARparStrand="--outSAMstrandField intronMotif"

echo "STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand"
STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand


samtools sort -@ 18 -o ./OUT/Aligned.sorted.out.bam ./OUT/Aligned.out.bam
samtools index ./OUT/Aligned.sorted.out.bam ./OUT/Aligned.sorted.out.bam.bai


mkdir ./OUT/Signal

echo "STAR --runMode inputAlignmentsFromBAM   --inputBAMfile ./OUT/Aligned.sorted.out.bam --outWigType bedGraph --outWigStrand Unstranded --outFileNamePrefix ./OUT/Signal/ --outWigReferencesPrefix chr"
STAR --runMode inputAlignmentsFromBAM   --inputBAMfile ./OUT/Aligned.sorted.out.bam --outWigType bedGraph --outWigStrand Unstranded --outFileNamePrefix ./OUT/Signal/ --outWigReferencesPrefix chr

###### bigWig conversion commands
grep ^chr ${GENPATH}/STARIndex/chrNameLength.txt > chrNL.txt

for imult in Unique UniqueMultiple
do
	grep ^chr ./OUT/Signal/Signal.$imult.str1.out.bg|sort -k1,1 -k2,2n > sig.tmp
	bedGraphToBigWig sig.tmp  chrNL.txt ./OUT/Signal/Signal.$imult.unstranded.bw

done

### HTSeq (genome based)
htseq-count -n 18 -f bam -r pos --idattr=gene_id --additional-attr=gene_name \
        ./OUT/Aligned.sorted.out.bam \
	${GENPATH}/gencode.v29.annotation.gtf \
        > ./OUT/Aligned.count.htseq.txt


mv OUT ./RNA-Seq/${LID}

done
