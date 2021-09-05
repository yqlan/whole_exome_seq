#!/bin/sh

cd /leofs/mishl_group/liangl/rare_dis/

awk -v FS="\t" '{print $1 ":" $2 ":" $3 ":" $4 ":" $5}' RG804_805_806_807_808Fam_anno.hg19_multianno_diffGenoType.txt > RG804_805_806_807_808Fam_anno.hg19_snp.txt

awk -v FS="\t" '{print $1 ":" $2 ":" $3 ":" $4 ":" $5}' RG798-799-800Fam_anno.hg19_multianno_diffGenoType.txt > RG798-799-800Fam_anno.hg19_snp.txt

awk -v FS="\t" '{print $1 ":" $2 ":" $3 ":" $4 ":" $5}' zhenfyFam_anno.hg19_multianno_diffGenoType.txt > zhenfyFam_anno.hg19_snp.txt

awk -v FS="\t" '{print $1 ":" $2 ":" $3 ":" $4 ":" $5}' RG801-802-803Fam_anno.hg19_multianno_diffGenoType.txt > RG801-802-803Fam_anno.hg19_snp.txt
