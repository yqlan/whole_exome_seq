#!/bin/bash

#PBS -N WGC811
#PBS -q largeq
#PBS -e WGC811_script_err
#PBS -o WGC811_script_log
#PBS -l mem=100gb,walltime=100:00:00,nodes=1:ppn=4
#HSCHED -s hschedd

#############################################################################################################################
cd /leofs/mishl_group/liangl/rare_dis
#ref=/leofs/mishl_group/liangl/rare_dis/hg19_gatkorder

#bwa=/software/biosoft/software/bwa-0.7.10/bwa
#java=/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java
#gatk=/software/biosoft/software/gatk/GenomeAnalysisTK.jar
#picard_dir=/software/biosoft/software/picard-tools-1.86
############################################################################################################################
#echo -n "the total num is\t" >> num.txt
#samtools view WGC034798U_dedup_realign_recal.bam | wc -l >> num.txt
#echo -n "the WGC num is\t" >> num.txt
#samtools view WGC034798U_dedup_realign_recal.bam | grep -o ":WGC" | wc -l >> num.txt

samtools view -h WGC034799U_dedup_realign_recal.bam | awk '{for(i=1;i<=NF;i++){if($i~/SM:/){print $i}}}' | sort -u >>WGC799_SM_PGinf.txt
samtools view -h WGC034799U_dedup_realign_recal.bam | awk '{for(i=1;i<=NF;i++){if($i~/PG:/){print $i}}}' | sort -u >>WGC799_SM_PGinf.txt

samtools view -h ../ZFY/WGC032239U.sorted_dedup_realign_recal.bam | awk '{for(i=1;i<=NF;i++){if($i~/SM:/){print $i}}}' | sort -u >>WGC239_SM_PGinf.txt
samtools view -h ../ZFY/WGC032239U.sorted_dedup_realign_recal.bam | awk '{for(i=1;i<=NF;i++){if($i~/PG:/){print $i}}}' | sort -u >>WGC239_SM_PGinf.txt

