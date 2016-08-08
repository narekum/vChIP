#!/bin/bash
script=$0
path=${script%/*}

name=$1
basename=${name%\.*}
dgf=$2
before=$3
after=$4
slice=$5
echo $basename
grep -P "\t\+\t" $name > $basename.fwd.bed
grep -P "\t\-\t" $name > $basename.rev.bed

$path/metagene_profiler_oriented.pl --exp $dgf --summ $basename.fwd.bed --out $basename.fwd.bed --bin 1 -bflank $before -aflank $after -ext 1 -print m -score z -rscore y -dcol 5 -th 5 --slice $slice
$path/metagene_profiler_oriented.pl --exp $dgf --summ $basename.rev.bed --out $basename.rev.bed --bin 1 -bflank $before -aflank $after -ext 1 -print m -score z -rscore y -dcol 5 -th 5 --slice $slice

perl -ne '$sum=0;chomp;@s=split(/\t/);@mat=@s[2..$#s];foreach (@mat) {$sum+=$_};$line=join "\t",$s[0],$s[1],$s[1]+1,".",,$sum,"+",@mat; print "$line\n" ' $basename.fwd.bed.fwd.zmatrix > $basename.fwd.zmatrix
perl -ne '$sum=0;chomp;@s=split(/\t/);@mat=@s[2..$#s];foreach (@mat) {$sum+=$_};$line=join "\t",$s[0],$s[1],$s[1]+1,".",,$sum,"+",@mat; print "$line\n" ' $basename.fwd.bed.rev.zmatrix > $basename.rev.zmatrix
perl -ne '$sum=0;chomp;@s=split(/\t/);@mat=@s[2..$#s];foreach (@mat) {$sum+=$_};$line=join "\t",$s[0],$s[1],$s[1]+1,".",,$sum,"+",@mat; print "$line\n" ' $basename.fwd.bed.both.zmatrix > $basename.both.zmatrix

perl -ne '$sum=0;chomp;@s=split(/\t/);@mat=@s[2..$#s];foreach (@mat) {$sum+=$_};$line=join "\t",$s[0],$s[1],$s[1]+1,".",,$sum,"-",@mat; print "$line\n" ' $basename.rev.bed.fwd.zmatrix >> $basename.fwd.zmatrix
perl -ne '$sum=0;chomp;@s=split(/\t/);@mat=@s[2..$#s];foreach (@mat) {$sum+=$_};$line=join "\t",$s[0],$s[1],$s[1]+1,".",,$sum,"-",@mat; print "$line\n" ' $basename.rev.bed.rev.zmatrix >> $basename.rev.zmatrix 
perl -ne '$sum=0;chomp;@s=split(/\t/);@mat=@s[2..$#s];foreach (@mat) {$sum+=$_};$line=join "\t",$s[0],$s[1],$s[1]+1,".",,$sum,"-",@mat; print "$line\n" ' $basename.rev.bed.both.zmatrix >> $basename.both.zmatrix

rm -f $basename.fwd.bed $basename.rev.bed $basename.fwd.bed.fwd.zmatrix $basename.fwd.bed.rev.zmatrix $basename.fwd.bed.both.zmatrix $basename.rev.bed.fwd.zmatrix $basename.rev.bed.rev.zmatrix $basename.rev.bed.both.zmatrix 

Rscript $path/scores.R $basename.fwd.zmatrix $basename.rev.zmatrix $basename.both.zmatrix $before $after 
