#!/usr/bin/perl 

#($slice,$column,$ascending,$master,$before,$after)=@ARGV;
($slice,$master,$percentage,$column,$ascending,$before,$after)=@ARGV;

$master=~/(.*).bed/;
$master=$1;

open BOTH,"$master.both.zmatrix";
chomp(@master_both=<BOTH>);
open FWD,"$master.fwd.zmatrix";
chomp(@master_fwd=<FWD>);
open REV,"$master.rev.zmatrix";
chomp(@master_rev=<REV>);

#print "$_\n" foreach @master_fwd ;

$slice=~/(.*).bed/;
$writedir=$1;
$writedir=join "_" , $writedir , "pc$percentage" , "col$column" , $ascending , "from-$before" , "to$after" ;
print "$writedir\n";
if ( -d $writedir ) {
	print "Directory $writedir already exists! exiting now\n";
	exit;
} else {
	mkdir $writedir ;
}

open OPEN,$slice;
foreach (<OPEN>){
	chomp;
	next if $_=~/^#/;
	push @slice,$_;
}
#chomp(@slice=<OPEN>);

if ($ascending eq "yes"){
	@slice_sorted= sort { (split /\t/,$a)[$column] <=> (split /\t/,$b)[$column]  } @slice ;
} else {
	@slice_sorted= sort { (split /\t/,$b)[$column] <=> (split /\t/,$a)[$column]  } @slice ;
}

#print "$_\n" foreach @slice ;
$total=$#slice_sorted + 1 ;
#print "$_\n" foreach @slice_sorted ;

#$percentage=20;
#$part=4;
#$writedir="./new";

for ($part=1;$part<=(100/$percentage);$part+=1) {
	my ($initial,$final)=returnPcSlice($total,$percentage,$part);
#	print "$total\t$initial\t$final\n";
	my @part_slice=@slice_sorted[$initial .. $final];

#print "$_\n" foreach @part_slice ;

	my @both=SLICE(\@master_both,\@part_slice);
	my @fwd=SLICE(\@master_fwd,\@part_slice);
	my @rev=SLICE(\@master_rev,\@part_slice);
	#print "$_\n" foreach @fwd ;


	open WR,">$writedir/$part.motif.txt";
	open REVWR,">$writedir/$part.revmotif.txt";
	foreach (@part_slice){
		my $motif=(split /\t/,$_)[6];
		my $revmotif = $motif ;
		$revmotif =~ tr/ACGT/TGCA/ ;
		$revmotif = reverse $revmotif ;
		print WR ">motif\n",$motif,"\n";
		print REVWR ">motif\n",$revmotif,"\n";
	}

	system "~/local/src/weblogo-3.3/weblogo -F png -o $writedir/$part.motif.png  < $writedir/$part.motif.txt";
	system "~/local/src/weblogo-3.3/weblogo -F png -o $writedir/$part.revmotif.png  < $writedir/$part.revmotif.txt";
	my ($fwdZscore,$posFwdPosMean,$posFwdCoordData)=stripCoordInfo(\@fwd);
	my ($revZscore,$posRevPosMean,$posRevCoordData)=stripCoordInfo(\@rev);
	my ($bothZscore,$posBothPosMean,$posRevCoordData)=stripCoordInfo(\@both);

	#print "$_\n" foreach @$posFwdPosMean ;
	$i=-$before;$j=$after;
	foreach ( @$posFwdPosMean ) { 
		$fwd{$i}{$part} = $_ ;
		$i++;
	}
	$i=-$before;
	foreach (@$posRevPosMean) {
		$rev{$i}{$part} = $_ ;
		$i++;
	}
	$i=-$before;
        foreach ( @$fwdZscore ) {
                $zfwd{$i}{$part} = $_ ;
                $i++;
        }
        $i=-$before;
        foreach (@$revZscore) {
                $zrev{$i}{$part} = $_ ;
                $i++;
        }
	$i=-$before;
	foreach (@$posBothPosMean) {
                $both{$i}{$part} = $_ ;
                $i++;
        }
	$i=-$before;
        foreach ( @$bothZscore ) {
                $zboth{$i}{$part} = $_ ;
                $i++;
        }
}

open FWD_WRITE,">$writedir/signature.fwd.dat";
print FWD_WRITE "#fwdposition\t";
print FWD_WRITE "$_\t" foreach ( 1 .. (100/$percentage) ) ;
print FWD_WRITE "\n";
foreach $pos ( -$before .. $after ) {
	print FWD_WRITE "$pos\t";
	foreach $div ( 1 .. (100/$percentage) ) {
		print FWD_WRITE $fwd{$pos}{$div} , "\t";
	}
	print FWD_WRITE "\n";
}

open REV_WRITE,">$writedir/signature.rev.dat";
print REV_WRITE "#revposition\t";
print REV_WRITE "$_\t" foreach ( 1 .. (100/$percentage) ) ;
print REV_WRITE "\n";
foreach $pos ( -$before .. $after ) {
        print REV_WRITE "$pos\t";
        foreach $div ( 1 .. (100/$percentage) ) {
                print REV_WRITE $rev{$pos}{$div} , "\t";
        }
        print REV_WRITE "\n";
}

open ZFWD_WRITE,">$writedir/zsignature.fwd.dat";
print ZFWD_WRITE "#zscorefwdposition\t";
print ZFWD_WRITE "$_\t" foreach ( 1 .. (100/$percentage) ) ;
print ZFWD_WRITE "\n";
foreach $pos ( -$before .. $after ) {
        print ZFWD_WRITE "$pos\t";
        foreach $div ( 1 .. (100/$percentage) ) {
                print ZFWD_WRITE $zfwd{$pos}{$div} , "\t";
        }
        print ZFWD_WRITE "\n";
}

open ZREV_WRITE,">$writedir/zsignature.rev.dat";
print ZREV_WRITE "#zscorerevposition\t";
print ZREV_WRITE "$_\t" foreach ( 1 .. (100/$percentage) ) ;
print ZREV_WRITE "\n";
foreach $pos ( -$before .. $after ) {
        print ZREV_WRITE "$pos\t";
        foreach $div ( 1 .. (100/$percentage) ) {
                print ZREV_WRITE $zrev{$pos}{$div} , "\t";
        }
        print ZREV_WRITE "\n";
}

open BOTH_WRITE,">$writedir/signature.both.dat";
print BOTH_WRITE "#bothposition\t";
print BOTH_WRITE "$_\t" foreach ( 1 .. (100/$percentage) ) ;
print BOTH_WRITE "\n";
foreach $pos ( -$before .. $after ) {
        print BOTH_WRITE "$pos\t";
        foreach $div ( 1 .. (100/$percentage) ) {
                print BOTH_WRITE $both{$pos}{$div} , "\t";
        }
        print BOTH_WRITE "\n";
}

open ZBOTH_WRITE,">$writedir/zsignature.both.dat";
print ZBOTH_WRITE "#zscorebothposition\t";
print ZBOTH_WRITE "$_\t" foreach ( 1 .. (100/$percentage) ) ;
print ZBOTH_WRITE "\n";
foreach $pos ( -$before .. $after ) { 
        print ZBOTH_WRITE "$pos\t";
        foreach $div ( 1 .. (100/$percentage) ) {
                print ZBOTH_WRITE $zboth{$pos}{$div} , "\t";
        }
        print ZBOTH_WRITE "\n";
}       


sub returnPcSlice {
        my ($total,$percentage,$part) =@_;
        my $percentage_number=int ( $percentage * ($total+1) / 100 );
        my $initial= ( $part - 1 ) * $percentage_number ;
        my $final= ( $part * $percentage_number ) - 1 ;
#       print "$percentage_number\t$part\t$initial\t$final\n";
        if ($final < $total ) {
                return ($initial,$final);
        } else {
                return (0,0);
        }
}

sub SLICE {
        my ($master,$slice)=@_;
        my %master;
        my %slice;
        my @return;

        #open OPEN,$master;
        foreach ( @$master ){
                chomp;
                next if $_=~/^#/;
                my @split=split (/\t/,$_);
                my $handle=join ":", $split[0] , $split[2] -1 , $split[2] , $split[5] ;
                $master{$handle}=$_;
                #print "$handle\n";
        }

        #open OPEN,$slice;
        foreach ( @$slice ){
                chomp;
                next if $_=~/^#/;
                my @split=split (/\t/,$_);
                my $handle;
                if ($split[5] eq "+"){
                        $handle=join ":", $split[0] , $split[1] , $split[1] + 1 , $split[5] ;
                } elsif ($split[5] eq "-"){
                        $handle=join ":", $split[0] , $split[2] -1 , $split[2] , $split[5] ;
                }
                push @return,$master{$handle} ;
                #print "$master{$handle}\n";
        }
        return @return ;
}

sub stripCoordInfo {
        my ($file)=@_;
        my $traindata;
        my $coordData;
	my @allcuts;
        my @trainfeatures;
        open FILE,$file;
        my @postionMean;
	my @postionZscore;
        my $n=0;
        foreach ( @$file){
                chomp; $n++;
                my @split=split /\t/;
                $traindata = join "\t", @split [ 6 .. $#split ] ;
                $coordData = join "\t", @split [ 0 .. 5 ] ;
                my $i=0;
                foreach ( @split [ 6 .. $#split ]  ){
                        $postionMean[$i] += $_ ;
			push @allcuts, $_ ;
                        $i++;
                }
                push @trainfeatures,$traindata;
                push @coordData,$coordData ;
        }
	my $mean=MEAN(\@allcuts);
	my ($stdev,$sterr)=STDEV(\@allcuts);
        for (my $i=0;$i<=$#postionMean;$i++){
		my $zscore= $postionMean[$i] /$n ;
                $postionMean[$i] = sprintf ("%.2f", $postionMean[$i] /$n ) ;
                $postionZscore[$i] = sprintf ("%.2f", ZSCORE($zscore,$stdev,$mean) ) ;
        }

        return (\@postionZscore,\@postionMean,\@coordData);
}

sub MEAN{
        my($data) = @_;
        if (not @$data) {
                die("Empty array $data\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = sprintf ("%.2f", $total / @$data ) ;
        return $average;
}
sub STDEV{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &MEAN($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = sprintf ( "%.2f", ($sqtotal / (@$data-1)) ** 0.5 ) ;
        my $err = sprintf ( "%.2f", ($std / sqrt(@$data)));
        return ($std,$err);
}

sub ZSCORE{
        my ($data,$stdev,$average)= @_;
        if ($stdev == 0 ){
                return (0);
        }
        my $zscore = ( $data - $average) / $stdev ;
        $zscore = sprintf ("%.2f", $zscore) ;
        return $zscore ;
}

