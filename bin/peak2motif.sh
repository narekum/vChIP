#! /bin/bash
export LC_ALL=C

script=$0
path=${script%/*}
pwd=`pwd`

#inputs
peakfile=$1
motiffile=$2
outdirname=$3
dgfcuts=$4

before=10
after=20
thres=1e-2
genomefasta=/mnt/WDMyBook2TB/library/indexes/genomes/bowtie2-indexes/hg19.fa
expchip=/mnt/WDMyBook2TB/NextGenSeq/human/alignments/K562.VEZF1.na.hg19.76bp.bt2.sensitive/bed/K562.VEZF1.na.hg19.76bp.bt2.sensitive.bam.uniq.bed
controlchip=/mnt/WDMyBook2TB/NextGenSeq/human/alignments/k562.IgG.na.hg19.76bp.bt2.sensitive/bed/k562.IgG.na.hg19.76bp.bt2.sensitive.bam.uniq.bed
#dgfcuts=/home/naren/sch/motif/fimo/coverage/coverage/positive_dataset/new/vChip/data/K562.DNaseDGF_SRR089677.normalized500bp.counts.VEZF1_76_bt2_sensitive_peaks.bed
chromHMM=/home/naren/sch/motif/fimo/coverage/coverage/positive_dataset/new/vChip/data/wgEncodeBroadHmmK562HMM.bed
gencode=/home/naren/sch/motif/fimo/coverage/coverage/positive_dataset/new/vChip/data/gencode.v19.annotation.gtf.tss.bed

if [ -d $pwd/$outdirname ]
then 
	echo "Directory $pwd/$outdirname already eitsts - Not doing anything further";
	exit;
else 
	mkdir $pwd/$outdirname
fi

echo Ranking peaks in $peakfile
sort -k4,4nr $peakfile | awk '{i++;print $1"\t"$2"\t"$3"\tRank_"i"\t"$4"\t+\t"$5}' | sort -k1,1 -k2,2n -k3,3n > $pwd/$outdirname/peak_ranked.bed
awk '{print $1"\t"$7"\t"$7+1"\t"$4"\t"$5"\t"$6}' $pwd/$outdirname/peak_ranked.bed > $pwd/$outdirname/peak_summits_ranked.bed

echo Getting fasta sequencing for peaks in $peakfile
bedtools getfasta -fi $genomefasta -bed $pwd/$outdirname/peak_ranked.bed -fo $pwd/$outdirname/peak_ranked.fasta

echo Searching for motif in $motiffile for instances in $peakfile
fimo --bgfile motif-file --max-seq-length 1000000000 --max-stored-scores 1000000000 --no-qvalue -oc $pwd/$outdirname/fimo_out --text --thresh $thres $motiffile $pwd/$outdirname/peak_ranked.fasta | $path/fimo2bed.pl > $pwd/$outdirname/peak_ranked.fimo.txt

echo Calculating $expchip read counts in identified motifs
bedtools coverage -a $expchip -b $pwd/$outdirname/peak_ranked.fimo.txt -counts > $pwd/$outdirname/peak_ranked.fimo.expcounts.txt
echo Calculating $controlchip read counts in identified motifs
bedtools coverage -a $controlchip -b  $pwd/$outdirname/peak_ranked.fimo.expcounts.txt -counts | sort -k1,1 -k2,2n -k3,3n > $pwd/$outdirname/peak_ranked.fimo.expcontrolcounts.txt
perl -e '$p="#chr\tstart\tend\tpeak_rank\tfimo-E-value\tstrand\tmotif\tname\tcorrected_counts\texp_counts\tcontrol_counts\tmacs_counts\tsummit";print "$p\n"'  > $pwd/$outdirname/peak_ranked.motif_master.bed
echo Intersection peak_ranked.bed with detected motifs
bedtools intersect -a $pwd/$outdirname/peak_ranked.fimo.expcontrolcounts.txt -b $pwd/$outdirname/peak_ranked.bed -wao | perl -ne '@sp=split;$sp[19]=$sp[8] - $sp[9];$p=join "\t", @sp[0..2,13,6,5,7,3,19,8,9,14,16];print "$p\n"' >> $pwd/$outdirname/peak_ranked.motif_master.bed
############
echo Calculating chromHMM state nearest to peak summit of $peakfile
bedtools intersect -a $pwd/$outdirname/peak_summits_ranked.bed -b $chromHMM -wao > $pwd/$outdirname/peak_summits_ranked_chromHMM.bed
echo Calculating  TSS distance nearest to peak summit of $peakfile
bedtools closest -a $pwd/$outdirname/peak_summits_ranked_chromHMM.bed -b $gencode -D a > $pwd/$outdirname/peak_summits_ranked_chromHMM_gencode.bed
rm $pwd/$outdirname/peak_summits_ranked_chromHMM.bed
echo combining chromHMM and nearest genCODE TSS with the motifs
$path/add_chromHMM_gencode.pl $pwd/$outdirname/peak_summits_ranked_chromHMM_gencode.bed $pwd/$outdirname/peak_ranked.motif_master.bed > $pwd/$outdirname/peak_ranked.motif_master_tss_gencode.bed
############
echo Calculating DGF signature for motifs
$path/rankSummits_unsorted_real.pl $pwd/$outdirname/peak_ranked.motif_master_tss_gencode.bed $dgfcuts > $pwd/$outdirname/motif_master.bed
$path/remove_motif_overlap.pl $pwd/$outdirname/motif_master.bed 22 MORE > $pwd/$outdirname/master_motifs_NR_DGF.bed
$path/remove_motif_overlap.pl $pwd/$outdirname/motif_master.bed 4 LESS > $pwd/$outdirname/master_motifs_NR_eValue.bed

rm -f $pwd/$outdirname/peak_ranked.fasta $pwd/$outdirname/peak_ranked.fimo.expcontrolcounts.txt $pwd/$outdirname/peak_ranked.fimo.expcounts.txt $pwd/$outdirname/peak_ranked.fimo.txt $pwd/$outdirname/peak_ranked.motif_master.bed $pwd/$outdirname/peak_ranked.motif_master_tss_gencode.bed $pwd/$outdirname/peak_summits_ranked_chromHMM_gencode.bed

echo  Calculating DGF cut matrix around the motif : before $before, after $after
$path/calc_distribution.sh $pwd/$outdirname/motif_master.bed $dgfcuts $before $after n 
$path/add_ranks.pl $pwd/$outdirname/motif_master.bed $pwd/$outdirname/motif_master.bed
$path/slice_matrix.pl $pwd/$outdirname/motif_master.bed $pwd/$outdirname/master_motifs_NR_DGF.bed
$path/slice_matrix.pl $pwd/$outdirname/motif_master.bed $pwd/$outdirname/master_motifs_NR_eValue.bed

echo Calculating scores at positions from -$before to $after and making graphs
Rscript $path/scores.R $pwd/$outdirname/master_motifs_NR_DGF.fwd.zmatrix $pwd/$outdirname/master_motifs_NR_DGF.rev.zmatrix $pwd/$outdirname/master_motifs_NR_DGF.both.zmatrix $before $after
Rscript $path/scores.R $pwd/$outdirname/master_motifs_NR_eValue.fwd.zmatrix $pwd/$outdirname/master_motifs_NR_eValue.rev.zmatrix $pwd/$outdirname/master_motifs_NR_eValue.both.zmatrix $before $after
Rscript $path/scores.R $pwd/$outdirname/motif_master.fwd.zmatrix $pwd/$outdirname/motif_master.rev.zmatrix $pwd/$outdirname/motif_master.both.zmatrix $before $after

echo Combining features for SVM 
$path/creat_training_dataset.pl $pwd/$outdirname/motif_master.bed > $pwd/$outdirname/motif_master_features.bed
$path/creat_training_dataset.pl $pwd/$outdirname/master_motifs_NR_DGF.bed > $pwd/$outdirname/master_motifs_NR_DGF_features.bed
$path/creat_training_dataset.pl $pwd/$outdirname/master_motifs_NR_eValue.bed > $pwd/$outdirname/master_motifs_NR_eValue_features.bed

exit

