#!/bin/bash

#PBS -N WGC034810U
#PBS -q largeq
#PBS -e WG034810U_script_err
#PBS -o WG034810U_script_log
#PBS -l mem=120gb,walltime=100:00:00,nodes=1:ppn=12
#HSCHED -s hschedd

#############################################################################################################################
cd /leofs/mishl_group/liangl/rare_dis
name=WGC034810U
R1=./Sample_${name}/${name}_combined_R1
R2=./Sample_${name}/${name}_combined_R2
ref=hg19_gatkorder
########################build index use bwa & samtools & picard#############################################################
#/software/biosoft/software/bwa-0.7.10/bwa index hg19_gatkorder.fa
#samtools faidx hg19_gatkorder.fa
#/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -jar /software/biosoft/software/picard-tools-1.86/CreateSequenceDictionary.jar R= hg19_gatkorder.fa O= hg19_gatkorder.dict

########################reads mapping use bwa###############################################################################
#/software/biosoft/software/bwa-0.7.10/bwa aln ${ref}.fa ${R1}.fastq > ${name}_R1.sai 2>>${name}_err

#/software/biosoft/software/bwa-0.7.10/bwa aln ${ref}.fa ${R2}.fastq > ${name}_R2.sai 2>>${name}_err

#/software/biosoft/software/bwa-0.7.10/bwa sampe -r "@RG\tID:WGC\tLB:WGC\tPL:illumina\tSM:WGC" ${ref}.fa ${name}_R1.sai ${name}_R2.sai ${R1}.fastq ${R2}.fastq > ${name}.sam 2>>${name}_err

#/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -jar /software/biosoft/software/picard-tools-1.86/SortSam.jar INPUT=${name}.sam OUTPUT=${name}.bam SORT_ORDER=coordinate TMP_DIR=${name}_tmp VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=false VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 CREATE_MD5_FILE=false 2>>${name}_err

#samtools index ${name}.bam 2>>${name}_err

#rm ${name}.sam

#######################Mark duplication use picard tools####################################################################
#/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -jar /software/biosoft/software/picard-tools-1.86/MarkDuplicates.jar INPUT=${name}.bam OUTPUT=${name}_dedup.bam METRICS_FILE=${name}_dedup.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=1 TMP_DIR=${name}_tmp VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=false PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 READ_NAME_REGEX=[a-zA-Z0-9]+:[0-9]:[0-9]+:[0-9]+:[0-9]+.* VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 CREATE_MD5_FILE=false 2>>${name}_err

#samtools index ${name}_dedup.bam 2>>${name}_err

#######################Indel Realignment use GATK###########################################################################
#/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${ref}.fa -I ${name}_dedup.bam -known: ./hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -known: ./hg19/1000G_phase1.indels.hg19.vcf -o ${name}_dedup.intervals 2>>${name}_err

#/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/gatk/GenomeAnalysisTK.jar -T IndelRealigner -R ${ref}.fa -I ${name}_dedup.bam -targetIntervals ${name}_dedup.intervals -known ./hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -known ./hg19/1000G_phase1.indels.hg19.vcf -o ${name}_dedup_realign.bam 2>>${name}_err

#samtools index ${name}_dedup_realign.bam 2>>${name}_err

#######################Base Recalibrator use GATK and picard################################################################
#/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/picard-tools-1.86/BamIndexStats.jar VALIDATION_STRINGENCY=SILENT INPUT=${name}.bam > ${name}_AlignPF.txt 2>>${name}_err

#/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/gatk/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${ref}.fa -I ${name}_dedup_realign.bam -knownSites: ./hg19/dbsnp_138.hg19.vcf -knownSites: ./hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -o ${name}_dedup_realign_recal.table 2>>${name}_err

/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/gatk/GenomeAnalysisTK.jar -T PrintReads -R ${ref}.fa -I ${name}_dedup_realign.bam -BQSR ${name}_dedup_realign_recal.table -o ${name}_dedup_realign_recal.bam 2>>${name}_err

samtools index ${name}_dedup_realign_recal.bam 2>>${name}_err


