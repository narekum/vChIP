#!/usr/bin/perl 
($file,$value_column,$which)=@ARGV;

#if $which eq "MORE" greater value motif is retained; 

open FILE,$file;
chomp(@allbedlines=<FILE>);
$header=shift @allbedlines ;

@allOrderedValidSummits=sort {(split /\t/,$a)[0] cmp (split /\t/,$b)[0] || (split /\t/,$a)[1] <=> (split /\t/,$b)[1] } @allbedlines ;
#print "$_\n" foreach @allOrderedValidSummits ;
#exit;
@survived=();
$first=shift @allOrderedValidSummits ;
#	print "first $first\n";
foreach $second  (@allOrderedValidSummits) {
#	print "second $second\n";
	@temp=();
	if ($which eq "MORE" ) {
		@temp=DECIDE_AND_RETURN_MORE($first,$second,$value_column);
	} elsif ($which eq "LESS") {
		@temp=DECIDE_AND_RETURN_LESS($first,$second,$value_column);
	}
#	print "each $_\n" foreach @temp;
	$first= pop @temp ;
	push @survived, @temp ;

}
print "$header\n";
print "$_\n" foreach @survived ;

sub DECIDE_AND_RETURN_MORE {
	my($first,$second,$value_column)=@_;
	my ($first_value, $first_start_coord,$first_end_coord, $second_value,$second_start_coord,$second_end_coord,@first,@second);
	@first=split(/\t/,$first);
	@second=split(/\t/,$second);

	if ( $first[0] ne $second[0] ) { # both motifs on separate chromosomes;
		return ($first,$second);
	}

	$first_value = $first[$value_column] ;
	$first_start_coord = $first[1] ;
	$first_end_coord = $first[2] ;

	$second_value = $second[$value_column] ;
	$second_start_coord = $second[1] ;
	$second_end_coord = $second[2] ;

#	print "$length, $first_value, $first_coord, $second_value, $second_coord  \n";

	if ( ( $second_start_coord - $first_end_coord ) >= 0 ) {
		return ( $first,$second);
	} else {

	        if ($first_value < $second_value ) { # if $first's evalue is less
	                return $second ; # $second is reported
	        } elsif ( $first_value > $second_value  ) { # if $seconds evlue is less
	                return $first ; # $first is reported
	        } else {
			return $first ; # return first is both are equal
	        }		
	}

}

sub DECIDE_AND_RETURN_LESS {
        my($first,$second,$value_column)=@_;
        my ($first_value, $first_start_coord,$first_end_coord, $second_value,$second_start_coord,$second_end_coord,@first,@second);
        @first=split(/\t/,$first);
        @second=split(/\t/,$second);

        if ( $first[0] ne $second[0] ) { # both motifs on separate chromosomes;
                return ($first,$second);
        }

        $first_value = $first[$value_column] ;
        $first_start_coord = $first[1] ;
        $first_end_coord = $first[2] ;

        $second_value = $second[$value_column] ;
        $second_start_coord = $second[1] ;
        $second_end_coord = $second[2] ;

#       print "$length, $first_value, $first_coord, $second_value, $second_coord  \n";

        if ( ( $second_start_coord - $first_end_coord ) >= 0 ) {
                return ( $first,$second);
        } else {

                if ($first_value < $second_value ) { # if $first's evalue is less
                        return $first ; # $first is reported
                } elsif ( $first_value > $second_value  ) { # if $seconds evlue is less
                        return $second ; # $second is reported
                } else {
                        return $second ; # return second is both are equal
                }
        }

}






#print "$_\n" for @allOrderedValidSummits ;
