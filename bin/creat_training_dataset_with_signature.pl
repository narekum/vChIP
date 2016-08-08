#!/usr/bin/perl
($positive,$signature_profile,$part)=@ARGV;

$positive=~s/.bed$//;
#$negative=~s/.bed$//;

($posFwd,$posFwdPosMean,$posFwdCoordData)=stripCoordInfo("$positive.fwd.zmatrix"); # Positive fwd strand,  Mean Counts
($posRev,$posRevPosMean,$posRevCoordData)=stripCoordInfo("$positive.rev.zmatrix");               # Positive rev strand,  Mean Counts

###########################################################
#($signature_profile,$part)=@ARGV;
$signature_profile=~/(.*).(fwd|rev|both).dat/;
$signature_profile=$1;
#print "$signature_profile\n";
$sigFwd=getProfile("$signature_profile.fwd.dat",$part);
$sigRev=getProfile("$signature_profile.rev.dat",$part);
$sigBoth=getProfile("$signature_profile.both.dat",$part);
#print "$_\n" foreach @$sigFwd ;
sub getProfile {
        my ($sigProfileFile,$part)=@_;
        open OPEN,$sigProfileFile;
        my @sig;
        while (<OPEN>){
                chomp;
                next if $_ =~ /^#/;
                my @split=split(/\t/,$_);
                push @sig,$split[$part];
        }
        return (\@sig);
}
##########################################################
$posFwdPosMean=$sigFwd;
$posRevPosMean=$sigRev;


@posFwdPosMean=@$posFwdPosMean;
@posRevPosMean=@$posRevPosMean;
@posFwd = @$posFwd ;
@posRev = @$posRev ;
@posFwdCoordData=@$posFwdCoordData;
@posRevCoordData=@$posRevCoordData;

#print "$_\n" foreach @$posRevPosMean;

#@posComb=combineFeatures(\@posFwd, \@posRev);

#($negFwd,$negFwdPosMean,$negFwdCoordData)=stripCoordInfo("$negative.fwd.zmatrix");
#($negRev,$negRevPosMean,$negFwdCoordData)=stripCoordInfo("$negative.rev.zmatrix");
#@negFwdPosMean=@$negFwdPosMean;
#@negRevPosMean=@$negRevPosMean;
#@negFwd=@$negFwd;
#@negRev=@$negRev;
#@negFwdCoordData=@$negFwdCoordData;
#@negFwdCoordData=@$negFwdCoordData;

#for ($i=0;$i<$#posRevPosMean;$i++){
	#print $posRevPosMean[$i] , "\t" , $posFwdPosMean[$i] , "\n";
	#print $i+1 , "\t",$posRevPosMean[$i] , "\t" , $posFwdPosMean[$i] , "\n";
#}

#@negComb=combineFeatures(\@negFwd, \@negRev);

#printFeatures(\@posComb,1);
#printFeatures(\@negComb,-1);

#&printFeaturesBoth($negFwd,$negFwdPosMean,$negFwdCoordData, $negRev,$negRevPosMean,1);
&printFeaturesBoth($posFwd,$posFwdPosMean,$posFwdCoordData, $posRev,$posRevPosMean,1);

sub printFeaturesBoth {
	my ($Fwd,$FwdPosMean,$FwdCoordData,$Rev, $RevPosMean, $classifyState)=@_;
	#my $svmFwd=returnFeatures($Fwd);
	#my $svmRev=returnFeatures($Rev);
	#my $svmFwdRev=returnFeature(@$Rev,);
	my @Fwd=@$Fwd; my @Rev=@$Rev; my @FwdCoordData=@$FwdCoordData;
	my @fwdCov;
	my @fwdMean;
	my @revCov;
	my @revMean;
	#my $FwdMean=join "\t" , @$FwdPosMean ;
	#my $RevMean=join "\t" , @$RevPosMean ;
	foreach (@$Fwd) {
		my @split=split(/\t/,$_);
		my $cov=covariance($FwdPosMean,\@split);
		push @fwdCov,$cov;
		my $fwdMean=MEAN(\@split);
		push  @fwdMean,$fwdMean;
		#print "cuts:$_ \nmean:$FwdMean\ncov:$cov\n\n";	
		#print "$cov\n";	
	}

	foreach (@$Rev){
		my @split=split(/\t/,$_);
		my $cov=covariance($RevPosMean,\@split);
		push @revCov,$cov;
		my $revMean=MEAN(\@split);	
		push  @revMean,$revMean;
	}
	#my $covarFwd=covariance($FwdPosMean,$Fwd);
	#print "$_\n" foreach @$svmFwd ;
	for (my $i=0;$i<=@$Fwd;$i++){
		print $classifyState,"\t",returnFeature( split(/\t/,$Fwd[$i]) , split(/\t/,$Rev[$i]) ),"\t#\t",$FwdCoordData[$i],"\t",$fwdCov[$i],"\t",$fwdMean[$i],"\t",sprintf("%.4f", $fwdCov[$i] * $fwdMean[$i]),"\t",$revCov[$i],"\t",$revMean[$i],"\t",sprintf("%.4f",$revCov[$i] * $revMean[$i] ) ,"\t",sprintf ("%.4f", ($fwdCov[$i] + $revCov[$i]) / 2 ) ,"\n" ;
		#print 
	}
	
}

#print "$_\n" foreach @pos;

sub stripCoordInfo {
	my ($file)=@_;
	my $traindata;
	my $coordData;
	my @trainfeatures;
	open FILE,$file;
	my @postionMean;
	my $n=0;
	while (<FILE>){
		chomp; $n++;
		my @split=split /\t/;
		$traindata = join "\t", @split [ 6 .. $#split ] ;
		$coordData = join "\t", @split [ 0 .. 5 ] ;
		my $i=0;
		foreach ( @split [ 6 .. $#split ]  ){
			$postionMean[$i] += $_ ;
			$i++;
		}
		push @trainfeatures,$traindata;
		push @coordData,$coordData ;
	}
	for (my $i=0;$i<=$#postionMean;$i++){
		$postionMean[$i] = sprintf ("%.2f", $postionMean[$i] /$n ) ;
	}
		
	return (\@trainfeatures,\@postionMean,\@coordData);	
}


sub combineFeatures {
	my ($fwd,$rev)=@_;
	my @fwd=@$fwd ;
	my @rev=@$rev ;
	#my $combineFeatures;
	my @combineFeatures;
	if ($#fwd != $#rev) {
		print "the two dataset to combine are not same in length\t";
		return;
	}

	for ( 0 .. $#fwd) {
		#print "$_\n";
		my $combineFeatures = join "\t" , $fwd[$_] , $rev[$_] ;
		push @combineFeatures, $combineFeatures ;
	}
	return (@combineFeatures) ;
}
sub returnFeature {
	my @features=@_;
	foreach ( 0 .. $#features){
		my $featurenumber= $_ + 1 ;
		$features[$_] = "$featurenumber:$features[$_]";
	}
	my $featureVector=join " ",@features;
	return $featureVector ;
}

sub returnFeatures {
	my ($features)=@_;
	my @features=@$features ;
	my @featuresReturn;
	foreach (@features){
		#my $featureVector .= "$classification ";	
		my @split=split /\t/ ;
		foreach ( 0 .. $#split) {
			my $featurenumber= $_ + 1 ;
			$split[$_] = "$featurenumber:$split[$_]";
			#$featureVector .= $featurenumber . ":" . $split[$_] . " "; 
			#print "$split[$_] ";
		}
		#print "\n";
		$featureVector=join " ",@split;
		push @featuresReturn,$featureVector;
		#print "$featureVector\n";
	}
	return (\@featuresReturn);
}

sub covariance {
        my ($array1ref, $array2ref) = @_; # array1 is the reference ; array2 is the variable
        my @array1ref=@$array1ref;
        my @array2ref=@$array2ref;
        if ($#array2ref != $#array1ref) {
                print "the two dataset to combine are not same in length in covariance subroutine\t";
                return;
        }
        my $array1mean=MEAN($array1ref);
        my $array2mean=MEAN($array2ref);
	#print "means $array1mean $array2mean\n";
        my $sumSqProduct= 0 ;
        for (my $i=0;$i<@array2ref;$i++){
                my $sqArray1 = ( $array1ref[$i] - $array1mean )  ;
                my $sqArray2 = ( $array2ref[$i] - $array2mean )  ;
                my $sqProduct = $sqArray1 * $sqArray2 ;
                $sumSqProduct += $sqProduct ;
        }
        my $cov = sprintf ("%.4f", $sumSqProduct / @array2ref );
        return ($cov);

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
        my $average = sprintf ("%.4f", $total / @$data ) ;
        return $average;
}








