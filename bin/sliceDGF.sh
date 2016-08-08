#!/bin/bash
dgf=$1
peaks=$2
length=$3

sort -k1,1 -k2,2n -k3,3n $peaks | bedtools slop -g ~/mount/publicdata/hg19/chrmSizes.hg19 -b $length | bedtools merge  -d 10  > $$.summits.sort.merge.bed
bedtools intersect -a $dgf -b $$.summits.sort.merge.bed -wa > dgf.bed
rm $$.summits.sort.merge.bed

