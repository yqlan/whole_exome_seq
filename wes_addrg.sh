#!/bin/bash

#PBS -N d_indel_reali
#PBS -q middleq
#PBS -e d_script_err
#PBS -o d_script_log
#PBS -l mem=60gb,walltime=24:00:00,nodes=1:ppn=8
#HSCHED -s hschedd

#############################################################################################################################

cd /leofs/mishl_group/liangl/rare_dis

#samtools view -h WG034804U_rg.bam > WG034804U_rg.sam
#samtools faidx hg19_gatkorder.fa
#java -jar /software/biosoft/software/picard-tools-1.86/CreateSequenceDictionary.jar R= hg19_gatkorder.fa O= hg19_gatkorder.dict

########################Add reads Group use picard tools####################################################################
#/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -jar /software/biosoft/software/picard-tools-1.86/AddOrReplaceReadGroups.jar I=WG034804U.bam O=WG034804U_rg.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=ture TMP_DIR=yq_tmp RGPL=illumina RGLB=WGC RGPU=WGC RGSM=WGC 2>>cmdline_err

########################Reorder Bam file use picard tools###################################################################
#/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -jar /software/biosoft/software/picard-tools-1.86/ReorderSam.jar VALIDATION_STRINGENCY=LENIENT TMP_DIR=yq_tmp CREATE_INDEX=ture I=WG034804U_rg.bam O=WG034804U_rg_gatkorder.bam R=hg19_gatkorder.fa 2>>cmdline_err

#samtools index WG034804U_rg_gatkorder.bam

#######################Mark duplication use picard tools####################################################################
#/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -jar /software/biosoft/software/picard-tools-1.86/MarkDuplicates.jar INPUT=WG034804U_rg_gatkorder.bam OUTPUT=WG034804U_rg_gatkorder_dedup.bam METRICS_FILE=WG034804U_rg_gatkorder_dedup.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=1 TMP_DIR=yq_tmp VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=true PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 READ_NAME_REGEX=[a-zA-Z0-9]+:[0-9]:[0-9]+:[0-9]+:[0-9]+.* VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 CREATE_MD5_FILE=false 2>>cmdline_err

#######################Indel Realignment use GATK###########################################################################
#/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R hg19_gatkorder.fa -I WG034804U_rg_gatkorder_dedup.bam -known: ./hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -known: ./hg19/1000G_phase1.indels.hg19.vcf -o WG034804U_rg_gatkorder_dedup.intervals 2>>cmdline_err

#/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/gatk/GenomeAnalysisTK.jar -T IndelRealigner -R hg19_gatkorder.fa -I WG034804U_rg_gatkorder_dedup.bam -targetIntervals WG034804U_rg_gatkorder_dedup.intervals -known ./hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -known ./hg19/1000G_phase1.indels.hg19.vcf -o WG034804U_rg_gatkorder_dedup_realign.bam 2>>cmdline_err

#######################Base Recalibrator use GATK and picard################################################################
/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/picard-tools-1.86/BamIndexStats.jar VALIDATION_STRINGENCY=LENIENT INPUT=WG034804U_rg_gatkorder.bam > WG034804U_rg_gatkorder_AlignPF.txt 2>>cmdline_err

/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/gatk/GenomeAnalysisTK.jar -T BaseRecalibrator -R hg19_gatkorder.fa -I WG034804U_rg_gatkorder_dedup_realign.bam -knownSites: ./hg19/dbsnp_138.hg19.vcf -knownSites: ./hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -o WG034804U_rg_gatkorder_dedup_realign_recal.table 2>>cmdline_err

/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/gatk/GenomeAnalysisTK.jar -T PrintReads -R hg19_gatkorder.fa -I WG034804U_rg_gatkorder_dedup_realign.bam -BQSR WG034804U_rg_gatkorder_dedup_realign_recal.table -o WG034804U_rg_gatkorder_dedup_realign_recal.bam 2>>cmdline_err

#######################Get Mutation use GATK################################################################################
#/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R hg19_gatkorder.fa -glm BOTH -I WGC_dedup_realign_recal.bam --dbsnp: dbsnp_138.hg19.vcf -stand_call_conf 30 -stand_emit_conf 10 -o WGC_dedup_realign_recal_rawU_snp_indel.vcf -minIndelCnt 4 -minIndelFrac 0.05 -dcov 1000

