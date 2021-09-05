#!/bin/bash

#PBS -N annovar798-799-800
#PBS -q middleq
#PBS -e annovar798-799-800_script_err
#PBS -o annovar798-799-800_script_log
#PBS -l mem=50gb,walltime=100:00:00,nodes=1:ppn=4
#HSCHED -s hschedd

#############################################################################################################################
cd /leofs/mishl_group/liangl/rare_dis

#perl /leofs/mishl_group/lanyq/software/annovar/table_annovar.pl /leofs/mishl_group/liangl/rare_dis/zhenfyFam_dedup_realign_recal_rawU_snpF_indelF.vcf /leofs/mishl_group/lanyq/work/rarediease/annovardb/ -buildver hg19 -out zhenfyFam_anno -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_sas,snp138,snp138Common,ljb26_all -operation g,r,r,f,f,f,f,f,f -nastring . -vcfinput 2>zhenfyFam_err

awk -v FS="\t" '/chr/{if($52=="PASS"){print $0}}' zhenfyFam_anno.hg19_multianno.txt > zhenfyFam_anno.hg19_multianno_pass.txt

awk -v FS="\t" '{if($6~/exonic|splicing/){print $0}}' zhenfyFam_anno.hg19_multianno_pass.txt > zhenfyFam_anno.hg19_multianno_func.refgene.txt

awk -v FS="\t" '{if($9!="."&&$9!="unknown"&&$9!="synonymous SNV"){print $0}}' zhenfyFam_anno.hg19_multianno_func.refgene.txt > zhenfyFam_anno.hg19_multianno_exonicfunc.refGene.txt

awk -v FS="\t" '{if($17=="."){print $0}}' zhenfyFam_anno.hg19_multianno_exonicfunc.refGene.txt > zhenfyFam_anno.hg19_multianno_snp138common.txt

awk -v FS="\t" -v OFS="\t" '{if($55~/\.\/\./||$56~/\.\/\./||$57~/\.\/\./){print $55,$56,$57,$0}}' zhenfyFam_anno.hg19_multianno_snp138common.txt > zhenfyFam_anno.hg19_multianno_noGenoType.txt

awk -v FS="\t" -v OFS="\t" '!/\.\/\./{split($55,case,":");split(case[2],casern,",");caserf=casern[2]/case[3];split($56,casef,":");split(casef[2],casefrn,",");casefrf=casefrn[2]/casef[3];split($57,casem,":");split(casem[2],casemrn,",");casemrf=casemrn[2]/casem[3];if(caserf>=0.6&&casefrf>=0.25&&casefrf<=0.6&&casemrf>=0.25&&casemrf<=0.6){print $0,case[1],caserf,case[3],casef[1],casefrf,casef[3],casem[1],casemrf,casem[3]}}' zhenfyFam_anno.hg19_multianno_snp138common.txt > zhenfyFam_anno.hg19_multianno_readref.txt

