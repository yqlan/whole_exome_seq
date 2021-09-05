#!/bin/perl

open OUT, '>/leofs/mishl_group/lim/ttttttt';
$file = `cat /leofs/mishl_group/lim/tmem_chr16snp.txt`;
print OUT "$file\n";
@num= $file=~m/HG00/g;
$num=@num;
print OUT "WGC number is $num\n";

