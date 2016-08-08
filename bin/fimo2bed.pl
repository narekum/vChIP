#!/usr/bin/perl 

#GATA1_6848	chr11:62166000-62166301	251	260	-	15.2241	5.14e-07	0.0388	AGAGATAAGA

#($fimo)=@ARGV;

#open OPEN,$fimo;
chomp(@fimo=<>);
#$out=FIMO2BED_center (\@fimo,"6" );
$out=FIMO2BED (\@fimo );
@out=@$out;

print "$_\n" foreach @out ;

sub FIMO2BED {
	my ($fimoformatedarray)=@_;
	my @fimoformatedarray=@$fimoformatedarray;
	my $each;
	my $output;
	my @output;
	foreach $each (@fimoformatedarray) {
		next if $each =~ /^#/ ;
		my @split=split(/\t/,$each) ;
		my $location=$split[1];
		$location=~/^(\S+):(\d+)-(\d+)/;
		my $chr=$1;
		my $start=$2;
		my $end=$3;
		my $motif_start=$start + $split[2] -1; # 0 based
		my $motif_end=$start + $split[3] ;       # 1 based
		$output= join("\t",$chr,$motif_start,$motif_end,join (",", $split[0],$split[6],$split[8]) ,$split[5],$split[4],@split[6 .. 8]);
		push @output,$output;
	}
	return (\@output);
}


sub FIMO2BED_center {
        my ($fimoformatedarray,$centeron)=@_;
        my @fimoformatedarray=@$fimoformatedarray;
        my $each;
        my $output;
        my @output;
        foreach $each (@fimoformatedarray) {
                my @split=split(/\t/,$each) ;
                my $location=$split[1];
                $location=~/^(\S+):(\d+)-(\d+)/;
                my $chr=$1;
                my $start=$2;
                my $end=$3;
		my $strand=$split[4];
		my $length=$split[3]-$split[2]+1;
		#my $centeron_rev=$length-$centeron;
		my $motif_start;
		#my $motif_end;
		if ($strand =~ /\+/ ){
                	$motif_start=$start + $split[2] -1 + $centeron -1 ; # 0 based
		} elsif ($strand =~ /-/){
			$motif_start=$start + $split[2] -1 + $length - $centeron;
		}
                my $motif_end=$motif_start+1 ;       # 1 based
		
		
                $output= join("\t",$chr,$motif_start,$motif_end,$split[0],$split[5],$split[4],@split[6 .. 8]);
                push @output,$output;
        }
        return (\@output);
}
