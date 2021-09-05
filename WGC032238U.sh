#!/bin/bash

#PBS -N zhenfyFam238
#PBS -q middleq
#PBS -e zhenfyFam238_script_err
#PBS -o zhenfyFam238_script_log
#PBS -l mem=60gb,walltime=100:00:00,nodes=1:ppn=15
#HSCHED -s hschedd

#############################################################################################################################
cd /leofs/mishl_group/liangl/rare_dis/
name=WGC032238U
R1=/leofs/mishl_group/liangl/rare_dis/Sample_${name}/${name}_combined_R1
R2=/leofs/mishl_group/liangl/rare_dis/Sample_${name}/${name}_combined_R2
ref=/leofs/mishl_group/liangl/rare_dis/hg19_gatkorder

bwa=/software/biosoft/software/bwa-0.7.10/bwa
java=/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java
gatk=/software/biosoft/software/gatk/GenomeAnalysisTK.jar
picard_dir=/software/biosoft/software/picard-tools-1.86
########################build index use bwa & samtools & picard#############################################################
#/software/biosoft/software/bwa-0.7.10/bwa index -a bwtsw hg19_gatkorder.fa
#samtools faidx hg19_gatkorder.fa
#/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -jar /software/biosoft/software/picard-tools-1.86/CreateSequenceDictionary.jar R= hg19_gatkorder.fa O= hg19_gatkorder.dict

########################reads mapping use bwa###############################################################################
#/software/biosoft/software/bwa-0.7.10/bwa aln ${ref}.fa ${R1}.fastq > ${name}_R1.sai 2>>${name}_err

/software/biosoft/software/bwa-0.7.10/bwa aln ${ref}.fa ${R2}.fastq > ${name}_R2.sai 2>>${name}_err

/software/biosoft/software/bwa-0.7.10/bwa sampe -r "@RG\tID:${name}\tLB:${name}\tPL:illumina\tSM:${name}" ${ref}.fa ${name}_R1.sai ${name}_R2.sai ${R1}.fastq ${R2}.fastq > ${name}.sam 2>>${name}_err

/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -jar /software/biosoft/software/picard-tools-1.86/SortSam.jar INPUT=${name}.sam OUTPUT=${name}.bam SORT_ORDER=coordinate TMP_DIR=${name}_tmp VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=false VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 CREATE_MD5_FILE=false 2>>${name}_err

samtools index ${name}.bam 2>>${name}_err

#rm ${name}.sam

#######################Mark duplication use picard tools####################################################################
/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -jar /software/biosoft/software/picard-tools-1.86/MarkDuplicates.jar INPUT=${name}.bam OUTPUT=${name}_dedup.bam METRICS_FILE=${name}_dedup.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=1 TMP_DIR=${name}_tmp VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=false PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 READ_NAME_REGEX=[a-zA-Z0-9]+:[0-9]:[0-9]+:[0-9]+:[0-9]+.* VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 CREATE_MD5_FILE=false 2>>${name}_err

samtools index ${name}_dedup.bam 2>>${name}_err

#######################Indel Realignment use GATK###########################################################################
/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${ref}.fa -I ${name}_dedup.bam -known: ./hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -known: ./hg19/1000G_phase1.indels.hg19.vcf -o ${name}_dedup.intervals 2>>${name}_err

/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/gatk/GenomeAnalysisTK.jar -T IndelRealigner -R ${ref}.fa -I ${name}_dedup.bam -targetIntervals ${name}_dedup.intervals -known ./hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -known ./hg19/1000G_phase1.indels.hg19.vcf -o ${name}_dedup_realign.bam 2>>${name}_err

samtools index ${name}_dedup_realign.bam 2>>${name}_err

#######################Base Recalibrator use GATK and picard################################################################
/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/picard-tools-1.86/BamIndexStats.jar VALIDATION_STRINGENCY=SILENT INPUT=${name}.bam > ${name}_AlignPF.txt 2>>${name}_err

/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/gatk/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${ref}.fa -I ${name}_dedup_realign.bam -knownSites: ./hg19/dbsnp_138.hg19.vcf -knownSites: ./hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -o ${name}_dedup_realign_recal.table 2>>${name}_err

/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java -Xmx45g -jar /software/biosoft/software/gatk/GenomeAnalysisTK.jar -T PrintReads -R ${ref}.fa -I ${name}_dedup_realign.bam -BQSR ${name}_dedup_realign_recal.table -o ${name}_dedup_realign_recal.bam 2>>${name}_err

samtools index ${name}_dedup_realign_recal.bam 2>>${name}_err

#######################Get Mutation use GATK################################################################################
#${java} -Xmx50g -jar ${gatk} -T UnifiedGenotyper -R ${ref}.fa -glm BOTH -I WGC032237U_dedup_realign_recal.bam -I WGC032238U_dedup_realign_recal.bam -I WGC032239U_dedup_realign_recal.bam --dbsnp: ./hg19/dbsnp_138.hg19.vcf -stand_call_conf 30 -stand_emit_conf 10 -o zhenfyFam_dedup_realign_recal_rawU_allvar.vcf -minIndelCnt 4 -minIndelFrac 0.05 -dcov 1000 2>>callvar_err

#${java} -Xmx50g -jar ${gatk} -T VariantRecalibrator -R ${ref}.fa -input zhenfyFam_dedup_realign_recal_rawU_allvar.vcf -mode SNP -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ./hg19/hapmap_3.3.hg19.vcf -resource:omni,known=false,training=true,truth=false,prior=12.0 ./hg19/1000G_omni2.5.hg19.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 ./hg19/1000G_phase1.indels.hg19.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./hg19/dbsnp_138.hg19.vcf -an DP -an FS -an MQRankSum -an ReadPosRankSum -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile zhenfyFam_dedup_realign_recal_rawU_snp.recal -tranchesFile zhenfyFam_dedup_realign_recal_rawU_snp.tranches -rscriptFile zhenfyFam_dedup_realign_recal_rawU_snp_plots.R 2>>callvar_err

#${java} -Xmx50g -jar ${gatk} -T VariantRecalibrator -R ${ref}.fa -input zhenfyFam_dedup_realign_recal_rawU_allvar.vcf -mode INDEL -resource:mills,known=true,training=true,truth=true,prior=12.0 ./hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -an DP -an FS -an MQRankSum -an ReadPosRankSum -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile zhenfyFam_dedup_realign_recal_rawU_indel.recal -tranchesFile zhenfyFam_dedup_realign_recal_rawU_indel.tranches -rscriptFile zhenfyFam_dedup_realign_recal_rawU_indel_plots.R 2>>callvar_err

#${java} -Xmx45g -jar ${gatk} -T ApplyRecalibration -R ${ref}.fa -input zhenfyFam_dedup_realign_recal_rawU_allvar.vcf -mode SNP --ts_filter_level 99.0 -recalFile zhenfyFam_dedup_realign_recal_rawU_snp.recal -tranchesFile zhenfyFam_dedup_realign_recal_rawU_snp.tranches -o zhenfyFam_dedup_realign_recal_rawU_snpF.vcf 2>>callvar_err

#${java} -Xmx45g -jar ${gatk} -T ApplyRecalibration -R ${ref}.fa -input zhenfyFam_dedup_realign_recal_rawU_snpF.vcf -mode INDEL --ts_filter_level 99.0 -recalFile zhenfyFam_dedup_realign_recal_rawU_indel.recal -tranchesFile zhenfyFam_dedup_realign_recal_rawU_indel.tranches -o zhenfyFam_dedup_realign_recal_rawU_snpF_indelF.vcf 2>>callvar_err

#######################annotate variants use annovar########################################################################
#perl /leofs/mishl_group/lanyq/software/annovar/table_annovar.pl /leofs/mishl_group/liangl/rare_dis/zhenfyFam_dedup_realign_recal_rawU_snpF_indelF.vcf /leofs/mishl_group/lanyq/work/rarediease/annovardb/ -buildver hg19 -out zhenfyFam_dedup_realign_recal_rawU_snpF_indelF_anno -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_sas,snp138,snp138Common,ljb26_all -operation g,r,r,f,f,f,f,f,f -nastring . -vcfinput 2>zhenfyFam_err

#awk -v FS="\t" '/chr/{if($52=="PASS"){print $0}}' zhenfyFam_dedup_realign_recal_rawU_snpF_indelF_anno.hg19_multianno.txt > zhenfyFam_dedup_realign_recal_rawU_snpF_indelF_anno.hg19_multianno_pass.txt

#awk -v FS="\t" '{if($6~/exonic|splicing/){print $0}}' zhenfyFam_dedup_realign_recal_rawU_snpF_indelF_anno.hg19_multianno_pass.txt > zhenfyFam_dedup_realign_recal_rawU_snpF_indelF_anno.hg19_multianno_func.refgene.txt

#awk -v FS="\t" '{if($9!="."&&$9!="unknown"&&$9!="synonymous SNV"){print $0}}' zhenfyFam_dedup_realign_recal_rawU_snpF_indelF_anno.hg19_multianno_func.refgene.txt > zhenfyFam_dedup_realign_recal_rawU_snpF_indelF_anno.hg19_multianno_exonicfunc.refGene.txt 

#awk -v FS="\t" '{if($17=="."){print $0}}' zhenfyFam_dedup_realign_recal_rawU_snpF_indelF_anno.hg19_multianno_exonicfunc.refGene.txt > zhenfyFam_dedup_realign_recal_rawU_snpF_indelF_anno.hg19_multianno_snp138common.txt

#awk -v FS="\t" -v OFS="\t" '{if($55~/\.\/\./||$56~/\.\/\./||$57~/\.\/\./){print $55,$56,$57,$0}}' zhenfyFam_dedup_realign_recal_rawU_snpF_indelF_anno.hg19_multianno_snp138common.txt > zhenfyFam_dedup_realign_recal_rawU_snpF_indelF_anno.hg19_multianno_noGenoType.txt

#awk -v FS="\t" -v OFS="\t" '!/\.\/\./{split($55,case,":");split(case[2],casern,",");caserf=casern[2]/case[3];split($56,casef,":");split(casef[2],casefrn,",");casefrf=casefrn[2]/casef[3];split($57,casem,":");split(casem[2],casemrn,",");casemrf=casemrn[2]/casem[3];if(caserf>=0.75&&casefrf>=0.25&&casefrf<=0.75&&casemrf>=0.25&&casemrf<=0.75){print $0,case[1],caserf,case[3],casef[1],casefrf,casef[3],casem[1],casemrf,casem[3]}}' zhenfyFam_dedup_realign_recal_rawU_snpF_indelF_anno.hg19_multianno_snp138common.txt > zhenfyFam_dedup_realign_recal_rawU_snpF_indelF_anno.hg19_multianno_readref.txt
