#!/usr/bin/perl
($master,$slice)=@ARGV;

$master=~/(.*).bed/;
$master=$1;

$slice=~/(.*).bed/;
$slice=$1;

@both=SLICE("$master.both.zmatrix","$slice.bed");
open BOTH, ">$slice.both.zmatrix";
print BOTH "$_\n" foreach @both ;

@fwd=SLICE("$master.fwd.zmatrix","$slice.bed");
open FWD, ">$slice.fwd.zmatrix";
print FWD "$_\n" foreach @fwd ;

@rev=SLICE("$master.rev.zmatrix","$slice.bed");
open REV, ">$slice.rev.zmatrix";
print REV "$_\n" foreach @rev ;


SLICE($master,$slice);

sub SLICE {
	my ($master,$slice)=@_;
	my %master;
	my %slice;
	my @return;

	open OPEN,$master;
	while (<OPEN>){
		chomp;
		next if $_=~/^#/;
		my @split=split (/\t/,$_);
		my $handle=join ":", $split[0] , $split[2] -1 , $split[2] , $split[5] ;
		$master{$handle}=$_;
		#print "$handle\n";
	}

	open OPEN,$slice;
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
		push @return,$master{$handle} ;
		#print "$master{$handle}\n";
	}
	return @return ;
}
