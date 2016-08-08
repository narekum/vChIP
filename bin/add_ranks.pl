#!/usr/bin/perl
($master,$slice)=@ARGV;

$master=~/(.*).bed/;
$master=$1;

$slice=~/(.*).bed/;
$slice=$1;

@both=SLICE("$master.bed","$slice.both.zmatrix");
open BOTH, ">$slice.both.zmatrix";
print BOTH "$_\n" foreach @both ;

@fwd=SLICE("$master.bed","$slice.fwd.zmatrix");
open FWD, ">$slice.fwd.zmatrix";
print FWD "$_\n" foreach @fwd ;

@rev=SLICE("$master.bed","$slice.rev.zmatrix");
open REV, ">$slice.rev.zmatrix";
print REV "$_\n" foreach @rev ;


SLICE($master,$slice);

sub SLICE {
	my ($master,$matrix)=@_;
	my %master;
	my %slice;
	my @return;

	open OPEN,$master;
	while (<OPEN>){
		chomp;
		next if $_=~/^#/;
		my @split=split (/\t/,$_);
		my $handle;
		if ($split[5] eq "+"){
			$handle=join ":", $split[0] , $split[1] , $split[1] + 1 , $split[5] ;
		} elsif ($split[5] eq "-"){
			$handle=join ":", $split[0] , $split[2] -1 , $split[2] , $split[5] ;
		}
		$master{$handle}=$split[3] ;
		#print "$master{$handle} \t$handle\n";
	}

        open OPEN,$matrix;
        while (<OPEN>){
                chomp;
                next if $_=~/^#/;
                my @split=split (/\t/,$_);
		$handle=join ":", $split[0] , $split[2] -1 , $split[2] , $split[5] ;
                #my $handle=join ":", $split[0] , $split[2] -1 , $split[2] , $split[5] ;
		$split[3] = $master{$handle} ;
		my $toreturn=join "\t" , @split ;
                push @return, $toreturn ;
                #print "$handle\n";
        }



	return @return ;
}
