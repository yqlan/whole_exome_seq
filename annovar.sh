perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
perl annotate_variation.pl -buildver hg19 -downdb genomicSuperDups humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2014oct humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp138 humandb/ 
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ljb26_all humandb/

perl /leofs/mishl_group/lanyq/software/annovar/table_annovar.pl /leofs/mishl_group/liangl/rare_dis/liukeweifam.filter.vcf /leofs/mishl_group/lanyq/work/rarediease/annovardb/ -buildver hg19 -out liukeweifam.filter.anno -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring NA -vcfinput
