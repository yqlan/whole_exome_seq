#!/bin/bash

#PBS -N alter_rg
#PBS -q middleq
#PBS -e alter_rg_script_err
#PBS -o alter_rg_script_log
#PBS -l mem=60gb,walltime=100:00:00,nodes=1:ppn=15
#HSCHED -s hschedd

#############################################################################################################################
cd /leofs/mishl_group/liangl/rare_dis
ref=/leofs/mishl_group/liangl/rare_dis/hg19_gatkorder

bwa=/software/biosoft/software/bwa-0.7.10/bwa
java=/software/biosoft/software/java1.7/jdk1.7.0_45/bin/java
gatk=/software/biosoft/software/gatk/GenomeAnalysisTK.jar
picard_dir=/software/biosoft/software/picard-tools-1.86

samtools view -h WGC034798U_dedup_realign_recal.bam | sed 's/:WGC/:WGC798/g' | samtools view -bS - > WGC798.bam
samtools index WGC798.bam
samtools view -h WGC034799U_dedup_realign_recal.bam | sed 's/:WGC/:WGC799/g' | samtools view -bS - > WGC799.bam
samtools index WGC799.bam
samtools view -h WGC034800U_dedup_realign_recal.bam | sed 's/:WGC/:WGC800/g' | samtools view -bS - > WGC800.bam
samtools index WGC800.bam

${java} -Xmx50g -jar ${gatk} -T UnifiedGenotyper -R ${ref}.fa -glm BOTH -I WGC798.bam -I WGC799.bam -I WGC800.bam --dbsnp: ./hg19/dbsnp_138.hg19.vcf -stand_call_conf 30 -stand_emit_conf 10 -o RG798-799-800Fam_dedup_realign_recal_rawU_allvar.vcf -minIndelCnt 4 -minIndelFrac 0.05 -dcov 1000 2>>RG798-799-800_err

${java} -Xmx50g -jar ${gatk} -T VariantRecalibrator -R ${ref}.fa -input RG798-799-800Fam_dedup_realign_recal_rawU_allvar.vcf -mode SNP -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ./hg19/hapmap_3.3.hg19.vcf -resource:omni,known=false,training=true,truth=false,prior=12.0 ./hg19/1000G_omni2.5.hg19.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 ./hg19/1000G_phase1.indels.hg19.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./hg19/dbsnp_138.hg19.vcf -an DP -an FS -an MQRankSum -an ReadPosRankSum -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile RG798-799-800Fam_dedup_realign_recal_rawU_snp.recal -tranchesFile RG798-799-800Fam_dedup_realign_recal_rawU_snp.tranches -rscriptFile RG798-799-800Fam_dedup_realign_recal_rawU_snp_plots.R 2>>RG798-799-800_err

${java} -Xmx50g -jar ${gatk} -T VariantRecalibrator -R ${ref}.fa -input RG798-799-800Fam_dedup_realign_recal_rawU_allvar.vcf -mode INDEL -resource:mills,known=true,training=true,truth=true,prior=12.0 ./hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -an DP -an FS -an MQRankSum -an ReadPosRankSum -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile RG798-799-800Fam_dedup_realign_recal_rawU_indel.recal -tranchesFile RG798-799-800Fam_dedup_realign_recal_rawU_indel.tranches -rscriptFile RG798-799-800Fam_dedup_realign_recal_rawU_indel_plots.R 2>>RG798-799-800_err

${java} -Xmx45g -jar ${gatk} -T ApplyRecalibration -R ${ref}.fa -input RG798-799-800Fam_dedup_realign_recal_rawU_allvar.vcf -mode SNP --ts_filter_level 99.0 -recalFile RG798-799-800Fam_dedup_realign_recal_rawU_snp.recal -tranchesFile RG798-799-800Fam_dedup_realign_recal_rawU_snp.tranches -o RG798-799-800Fam_dedup_realign_recal_rawU_snpF.vcf 2>>RG798-799-800_err

${java} -Xmx45g -jar ${gatk} -T ApplyRecalibration -R ${ref}.fa -input RG798-799-800Fam_dedup_realign_recal_rawU_snpF.vcf -mode INDEL --ts_filter_level 99.0 -recalFile RG798-799-800Fam_dedup_realign_recal_rawU_indel.recal -tranchesFile RG798-799-800Fam_dedup_realign_recal_rawU_indel.tranches -o RG798-799-800Fam_dedup_realign_recal_rawU_snpF_indelF.vcf 2>>RG798-799-800_err
