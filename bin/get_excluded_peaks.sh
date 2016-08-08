#!/bin/bash

DGFfull=$1
TFfull=$2
outdirname=$3
DGF=${DGFfull##*/}
TF=${TFfull##*/}
pwd=`pwd`

dgfcuts=/home/naren/sch/motif/fimo/coverage/coverage/positive_dataset/new/vChip/data/K562.DNaseDGF_SRR089677.normalized500bp.counts.VEZF1_76_bt2_sensitive_peaks.bed
expchip=/mnt/WDMyBook2TB/NextGenSeq/human/alignments/K562.VEZF1.na.hg19.76bp.bt2.sensitive/bed/K562.VEZF1.na.hg19.76bp.bt2.sensitive.bam.uniq.bed

if [ -d $pwd/$outdirname ]
then
        echo "Directory $pwd/$outdirname already eitsts - Not doing anything further";
        exit;
else
        mkdir $pwd/$outdirname
fi


#bedtools intersect -a $DGFfull -b $TFfull -loj -v > $DIRNAME/$DGF.wo.$TF.intersect.bed
#awk '{if ($11 == ".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10  }' $DGF.$TF.intersect.bed > $DGF.wo.$TF.bed
#./$DGF.$TF.intersect.bed
#sort -R $DGF.wo.$TF.bed | head -`cat $TF | wc -l` > $DGF.wo.$TF.bed

bedtools intersect -a $DGFfull -b $TFfull -loj -v  | awk '{h=($3-$2)/2;;printf "%s\t%.0f\t%.0f\t%s\t%.0f\n",  $1, $2, $3,".",$2+h}' > $pwd/$outdirname/$DGF.wo.$TF.intersect.bed
bedtools coverage -a $expchip -b $pwd/$outdirname/$DGF.wo.$TF.intersect.bed -counts | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$5}' | sort -k4n > $pwd/$outdirname/$DGF.wo.$TF.intersect_sort.bed
rm $pwd/$outdirname/$DGF.wo.$TF.intersect.bed 
