#!/bin/bash

#PBS -N annovar804_805_806_807_808
#PBS -q middleq
#PBS -e annovar804_805_806_807_808_script_err
#PBS -o annovar804_805_806_807_808_script_log
#PBS -l mem=50gb,walltime=100:00:00,nodes=1:ppn=4
#HSCHED -s hschedd

#############################################################################################################################
cd /leofs/mishl_group/liangl/rare_dis

#perl /leofs/mishl_group/lanyq/software/annovar/table_annovar.pl /leofs/mishl_group/liangl/rare_dis/RG804_805_806_807_808Fam_dedup_realign_recal_rawU_snpF_indelF.vcf /leofs/mishl_group/lanyq/work/rarediease/annovardb/ -buildver hg19 -out RG804_805_806_807_808Fam_anno -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_sas,snp138,snp138Common,ljb26_all -operation g,r,r,f,f,f,f,f,f -nastring . -vcfinput 2>RG804_805_806_807_808Fam_err

#awk -v FS="\t" '/^chr/{if($52=="PASS"){print $0}}' RG804_805_806_807_808Fam_anno.hg19_multianno.txt > RG804_805_806_807_808Fam_anno.hg19_multianno_pass.txt

#awk -v FS="\t" '{if($6~/exonic|splicing/){print $0}}' RG804_805_806_807_808Fam_anno.hg19_multianno_pass.txt > RG804_805_806_807_808Fam_anno.hg19_multianno_func.refgene.txt

#awk -v FS="\t" '{if($9!="."&&$9!="unknown"&&$9!="synonymous SNV"){print $0}}' RG804_805_806_807_808Fam_anno.hg19_multianno_func.refgene.txt > RG804_805_806_807_808Fam_anno.hg19_multianno_exonicfunc.refGene.txt

#awk -v FS="\t" '{if($17=="."){print $0}}' RG804_805_806_807_808Fam_anno.hg19_multianno_exonicfunc.refGene.txt > RG804_805_806_807_808Fam_anno.hg19_multianno_snp138common.txt

#awk -v FS="\t" -v OFS="\t" '{if($55~/\.\/\./||$56~/\.\/\./||$57~/\.\/\./||$58~/\.\/\./||$59~/\.\/\./){print $55,$56,$57,$58,$59,$0}}' RG804_805_806_807_808Fam_anno.hg19_multianno_snp138common.txt > RG804_805_806_807_808Fam_anno.hg19_multianno_noGenoType.txt

#awk -v FS="\t" -v OFS="\t" '!/\.\/\./{split($55,case,":");split(case[2],casern,",");caserf=casern[2]/case[3];split($56,casef,":");split(casef[2],casefrn,",");casefrf=casefrn[2]/casef[3];split($57,casem,":");split(casem[2],casemrn,",");casemrf=casemrn[2]/casem[3];split($58,cases,":");split(cases[2],casesrn,",");casesrf=casesrn[2]/cases[3];split($59,caseb,":");split(caseb[2],casebrn,",");casebrf=casebrn[2]/caseb[3];if(caserf>=0.6&&casefrf>=0.15&&casefrf<=0.6&&casemrf>=0.15&&casemrf<=0.6&&casesrf>=0.15&&casesrf<=0.6&&casebrf>=0.15&&casebrf<=0.6){print $0,case[1],caserf,case[3],casef[1],casefrf,casef[3],casem[1],casemrf,casem[3],cases[1],casesrf,cases[3],caseb[1],casebrf,caseb[3]}}' RG804_805_806_807_808Fam_anno.hg19_multianno_snp138common.txt > RG804_805_806_807_808Fam_anno.hg19_multianno_readref.txt

awk -v FS="\t" -v OFS="\t" '!/\.\/\./{split($55,case,":");split(case[2],casern,",");caserf=casern[2]/case[3];split($56,casef,":");split(casef[2],casefrn,",");casefrf=casefrn[2]/casef[3];split($57,casem,":");split(casem[2],casemrn,",");casemrf=casemrn[2]/casem[3];split($58,cases,":");split(cases[2],casesrn,",");casesrf=casesrn[2]/cases[3];split($59,caseb,":");split(caseb[2],casebrn,",");casebrf=casebrn[2]/caseb[3];if(case[1]!=casef[1]&&case[1]!=casem[1]&&case[1]!=cases[1]&&case[1]!=caseb[1]){print $0}}' RG804_805_806_807_808Fam_anno.hg19_multianno_snp138common.txt > RG804_805_806_807_808Fam_anno.hg19_multianno_readref.txt
