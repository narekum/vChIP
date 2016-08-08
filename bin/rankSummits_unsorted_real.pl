#! /usr/bin/perl
($summitFile,$expBedfile)=@ARGV;

$readType="exp";
open READS,$expBedfile;
while (<READS>){
        chomp;
	@split=split(/\t/,$_);
	#$$chrhash{$split[1]}+=$split[$dcol - 1] ;
        #$_=~/^(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\d+)\t(\S+)/;
	$chrhash=$readType.$split[0]  ;
        if ( $split[5] eq "-") {
                $$chrhash{ $split[2] - 1 }+=$split[4];
        } else {
                $$chrhash{$split[1]}+=$split[4];
        }
}

open PEAKS,$summitFile;
while (<PEAKS>){
	chomp;
	$line=$_;
	$line=~/^(\S+)\t(\d+)\t(\d+)\t(\S+)/;
	if ($line=~/^#/) { print "$line\tMotifDGF\tPreviousDGF\tNextDGF\tDGFIndex\n"; next;}
	$chr=$1;
	$start=$2;
	$end=$3;
	$length=$end-$start;
		$expCount=reportMean($chr,$start,$end,"exp");
		$expCount=sprintf ("%.2f",$expCount);
		$previous=reportMean($chr,$start - $length ,$start,"exp");
		$previous=sprintf ("%.2f",$previous);
		$next=reportMean($chr,$end ,$end + $length ,"exp");
		$next=sprintf ("%.2f",$next);
		#$dfg_index=sprintf("%.2f", (($previous + $next) / ( 2*$expCount )));
		$dfg_index=sprintf("%.2f", (( $previous + $next - 2*$expCount ) / 2 ));
		print "$line\t$expCount\t$previous\t$next\t$dfg_index\n";
		#push @unsortedList, join "\t", $chr,$summits,$summits + 1,$expCount,$controlCount,$correctedCount,$expRPM,$controlRPM,$correctedCountRPM ;
		#printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$chr,$summits,$summits + 1,$expCount,$controlCount,$correctedCount,$expRPM,$controlRPM,$correctedCountRPM ;
}

sub reportMean {
        my ($chr,$start,$end,$exp) = @_ ;
        my $coordinates ;
        my $count=0 ;
	$chr=$exp.$chr;
        for $coordinates ( $start .. ( $end - 1 )  ) {
                if (exists $$chr{$coordinates} ) {
                        $count += $$chr{$coordinates} ;
#                        print "$coordinates\texists=", $$chr{$coordinates} ,"\n" ;
                } else {
 #                       print "$coordinates\tnope\n" ;
                }
        }

        #print "total counts= $count \n";
        return $count ;

}

