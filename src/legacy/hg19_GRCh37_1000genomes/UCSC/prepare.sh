#!/usr/bin/env bash

creation of additional tracks

zcat CpGIslands_plain.bed.gz | perl -ne 'if($_ =~ /#/){print};chomp;@a=split;print "$a[0]\t",$a[1]-2000,"\t$a[1]\n";print "$a[0]\t$a[2]\t",$a[2]+2000, "\n"' | sort -k1,1V -k2,2n -k3,3n | bgzip > CpGIslands_shores2000bp_plain.bed.gz

zcat RefSeq_Sept_2013_from_annovar_Genes_plain.bed.gz | perl /home/hutter/Perlprog/ucscDB/PromoRegfromMerged.pl - 2000 500 | sort -k1,1V -k2,2n -k3,3n | bgzip > RefSeq_Sept_2013_from_annovar_TSS_2000up_500down_plain.bed.gz

