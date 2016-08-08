#! /usr/bin/perl 

#############################################
##  Narendra Kumar, PhD                     #
##  Dr Adam West Lab                        #
##  Epigenetics Unit                        #
##  Institute of Cancer Sciences            #
##  University of Glasgow, UK               #
##  narekum@gmail.com                       #
#############################################

#use strict ;
use Getopt::Long;
use LWP::Simple;
use Parallel::ForkManager;

my $binwidth=100;            # --binwidth     : size of each bin in bp (default 100)
my $bflank=1000;             # --bflank        : window on each side on the reference peak positions (default 1000)
my $aflank=1000;             # --aflank        : window on each side on the reference peak positions (default 1000)
my $extend_length=150;       # --extendlength : how much the reads should be extended (default 100)
my $control="NULL";          # --control      : control reads file (required)
my $nproc=2;                 # --thread       : no of processors (default 2)
my $toprint="s";             # --print        : what to report in the output s=scores; m=matrix; b=both (default s)
my $mean_or_zscore="z";      # --score        : score output type m=mean; z=zscore (default z)
my $rscore="n";              # --rscore       : if the read scores are to be used y=yes; n=no (default n)
my $dcol=5;                  # --dcol         : data column in the bed file to be used if --rscore=y ; (default 5)
my $strand="fwd";            # --strand       : which strand to use for counting aligned reads --strand=fwd; --strand=rev (default fwd)
my $slice_or_not="y";        # --slice        : wheather or not take into memory only reads in summits regions y=yes; n=no (default no)
my $normalize=0;             # --normalize    : 1000000 for reads per million (default no normalization)
my $groupsize=0;             # --groupsize    : no to summits to devide into groups (default no division)

my $exp;
my $summitfile;
my $outfile;
my $time_tag;
my $start_time;
my @bincontainer;
my $totalExpReads;
my $totalCtrReads;
my $chr;
my $groups;
my %chr;
my %binmatrix;
my @outmatrix;
my @outscores;
my $outmatrix;
my $outscores;
my $binzmatrix;
my $binscoresd;
my $chrhash;

$time_tag=$start_time=time;

GetOptions ("binwidth=i" => \$binwidth,
            "bflank=i" => \$bflank,
            "aflank=i" => \$aflank,
            "exp=s" => \$exp,
            "control=s" => \$control,
            "summitfile=s" => \$summitfile,
            "extendlength=i" => \$extend_length,
            "thread=i" => \$nproc,
            "outfile=s" => \$outfile,
            "print=s" => \$toprint,
            "scoretype=s" => \$mean_or_zscore,
            "rscore=s" => \$rscore,
            "dcol=i" => \$dcol,
            "strand=s" => \$strand,
            "slice=s" => \$slice_or_not,
            "normalize=i" => \$normalize,
	    "groupsize=i" => \$groupsize
           )
or die("Error in command line arguments\n");

if ( $exp eq "" || $summitfile eq "" || $outfile eq ""  ) {
	print "

ERROR: Required parameters missing: required parameters are --exp, --summitfile and --outfile

Usage :metagene_profiler_oriented.pl --exp <exp bed file>  --summitfile <summit bed file> --outfile <outfile>

Options:
	--binwidth     : size of each bin in bp (default 100)
	--bflank       : window before (upstream ) the reference peak positions (default 1000)
	--aflank       : window after (downstrem) the reference peak positions (default 1000)
	--extendlength : how much the reads should be extended (default 100)
	--control      : control reads file (required)
	--thread       : no of processors (default 2)
	--print        : what to report in the output s=scores; m=matrix; b=both (default s)
	--score        : score output type m=mean; z=zscore (default z)
	--rscore       : if the read scores are to be used y=yes; n=no (default n)
	--dcol         : data column in the bed file to be used if --rscore=y ; (default 5)
	--strand       : which strand to use for counting aligned reads --strand=fwd; --strand=rev (default fwd)
	--slice        : wheather or not take into memory only reads in summits regions y=yes; n=no (default no)
	--normalize    : 1000000 for reads per million (default no normalization)
	--groupsize    : no to summits to devide into groups (default no division)
";

	exit;
}

print "Calculating bins\n";
@bincontainer = BINS($binwidth,$bflank,$aflank); #PRINTTIME("Bin calculated");
#print "bin $_\n" foreach @bincontainer;
PRINTTIME ("Calculating bins");

#print "read $_\n" foreach @$slice;

if ($slice_or_not eq "y" ) {
	print "Slicing exp reads around summits\n";
	$slice=SLICE($exp,$summitfile,$bflank,$aflank);
	print "Reading aligned tags\n";
	PRINTTIME ("Slicing exp reads around summits");
	$totalExpReads=READALIGNEDSLICEDTAGS($slice,$exp,$extend_length,$rscore,$dcol);
	if ($control ne "NULL") {
		$totalCtrReads=READALIGNEDSLICEDTAGS($control,$extend_length,$rscore,$dcol);
	}
} elsif ($slice_or_not eq "n") {
	$totalExpReads=READALIGNEDTAGS($exp,$extend_length,$rscore,$dcol);
} else {
	print "unknown --slice $slice_or_not value\n";
}
PRINTTIME("Reading aligned tags");

print "Reading summits\n";
($chr,$groups)=READSUMMITS($summitfile,$groupsize); #PRINTTIME("SUMMITS READ");
%chr=%$chr;
PRINTTIME("Reading Summits");

&EXECUTE (\%chr,\@bincontainer,$exp,$control,$binwidth,$extend_length,$normalize,$totalExpReads,$totalCtrReads,$binwidth,"fwd",$toprint,$groups,$outfile, $mean_or_zscore);
&EXECUTE (\%chr,\@bincontainer,$exp,$control,$binwidth,$extend_length,$normalize,$totalExpReads,$totalCtrReads,$binwidth,"rev",$toprint,$groups,$outfile, $mean_or_zscore);
&EXECUTE (\%chr,\@bincontainer,$exp,$control,$binwidth,$extend_length,$normalize,$totalExpReads,$totalCtrReads,$binwidth,"both",$toprint,$groups,$outfile, $mean_or_zscore);
#%binmatrix=RETURNBINCOUNTSALLCHR(\%chr,\@bincontainer,$exp,$control,$binwidth,$extend_length,$normalize,$totalExpReads,$totalCtrReads,$binwidth);
sub EXECUTE {
	my ($chr,$bincontainer,$exp,$control,$binwidth,$extend_length,$normalize,$totalExpReads,$totalCtrReads,$binwidth,$strand,$toprint,$groups,$outfile,$mean_or_zscore)=@_;
	my %binmatrix=();
	my ($outscores,$outmatrix,$binzmatrix,$binscoresd,$outscores);
	print "Calculating distribution\n";
	%binmatrix=RETURNBINCOUNTSALLCHR($chr,$bincontainer,$exp,$control,$binwidth,$extend_length,$normalize,$totalExpReads,$totalCtrReads,$binwidth,$strand);
	PRINTTIME("Calculating distribution");

	my @outmatrix=(); # array to store matrix. To be printed later on
	my @outscores=(); # array to store scores. To be printed later on
	print "Printing results\n";
	if ($toprint eq "m" || $toprint eq "M" || $toprint eq "b" || $toprint eq "B") {
		$outmatrix=\@outmatrix;
	foreach  (@$groups) {
		#print "each groups ---- $_\n";
		if ( $groupsize != 0 ) {   # enter this condition only if data is splitted into groups, in that case we dont want double counts
			next if $_ =~ /all/ ; # Do not add the first element (all summits) as that would lead to double count
		}
		$outmatrix = PRINTMATRIX($chr,$bincontainer,\%binmatrix,$outmatrix,$_);
	}


	PRINTARRAY($outmatrix,"$outfile.$strand.zmatrix");

	} 

	if ( $toprint eq "s" || $toprint eq "S" || $toprint eq "b" || $toprint eq "B" ) {
		$outscores=\@outscores;
	#print "$toprint------$outfile.$strand.zscores\n";
		foreach  (@$groups) {
	#		my ($binzmatrix,$binscoresd)=BINZSOCRES(\%chr,\@bincontainer,\%binmatrix,$mean_or_zscore,$_);
			($binzmatrix,$binscoresd)=BINZSOCRES($chr,$bincontainer,\%binmatrix,$mean_or_zscore,$_);

			$outscores = REPORTRESULTS ($binzmatrix,$binscoresd,$bincontainer,$outscores,$_);
			PRINTARRAY($outscores,"$outfile.$strand.zscores");
		}
	}
PRINTTIME("Printing results");

}

#########################################################
sub SLICE {
        my ($reads,$summits,$before,$after) = @_ ;
        my (@summ,$summ,@expanded,$slice);
        open SUMM, $summits;
        chomp(@summ=<SUMM>);
        close SUMM;
        foreach $summ (@summ){
                my @sp = split (/\t/,$summ);
                my $expanded;
                if ($sp[5] eq "+" ){
                        push @expanded, join ( "\t", $sp[0] , $sp[1] - $before , $sp[2] + $after,".","0","+"  ) ;
#                       print "$expanded\n" ;
                } elsif ( $sp[5] eq "-" ) {
                        push @expanded, join ( "\t", $sp[0] , $sp[1] - $after , $sp[2] + $before,".","0","-"  ) ;
                }
        }
        open WRITE,">$$.summits";
        print WRITE "$_\n" foreach @expanded ;
	print "sorting expanded sliced summits\n";
        system ("sort -k1,1 -k2,2n -k3,3n $$.summits > $$.summits.sort.bed");
	PRINTTIME ("sorting expanded sliced summits");
	print "merging nerby summits for extracting reads form exp file\n";
        system ("bedtools merge -i $$.summits.sort.bed -d 10 > $$.summits.sort.merge.bed");
	print ("merging nerby summits for extracting reads form exp file");
	print "extracting reads from exp bedfile in expanded merged summit regions\n";
        chomp (@slice=`bedtools intersect -a $reads -b $$.summits.sort.merge.bed -wa`) ;
        system ("bedtools intersect -a $reads -b $$.summits.sort.merge.bed -wa > VEZF1.slice.300each_allmotifs_norm_counts.bed") ;
	print ("extracting reads from exp bedfile in expanded merged summit regions");
        system ("rm -f $$.summits.sort.bed $$.summits.sort.merge.bed $$.summits");
	print "MASTER PRINTED\n";
        return (\@slice);
	
}

sub PRINTARRAY {
	my ($array,$outfile)=@_;
	open WRITE, ">$outfile" ;
	foreach (@$array){
		print WRITE "$_\n";
	}
	close WRITE ;
}

sub PRINTMATRIX {
        my ($chr,$bincontainer,$binmatrix,$outmatrix,$group)=@_;
        my @bincontainer=@$bincontainer;
	#print "$_\t" foreach @bincontainer ; print "\n";
        my %binmatrix=%$binmatrix;
        my %chr = %$chr;
        my ($coor,$keys,$bin,$groupchr,$groupcoord);
	my $outline ;
	my @outmatrix=@$outmatrix;
	foreach ( @$group ) {
		($groupchr,$groupcoord) = split(/-/,$_);
		$outline =  "$groupchr\t$groupcoord" ;
		foreach $bin ( @bincontainer ) {
			#print "i--$bin\t";
			$outline = $outline . "\t" . $binmatrix{$groupchr}{$groupcoord}{$bin} ;
			#print "->", $binmatrix{$groupchr}{$groupcoord}{$bin} ,"<-" ;print "\n";
		}
		push @outmatrix, $outline ;
	}
	return (\@outmatrix);
}

sub REPORTRESULTS {
	my ($binzmatrix,$binscoresd,$bincontainer,$outarray,$group)=@_;
	my %binzmatrix = %$binzmatrix ;
	my %binscoresd = %$binscoresd ;
	my @bincontainer=@$bincontainer;
	my($bin,$i);
	my @outarray=@$outarray;
	if (scalar @outarray == 0 ) {
		push @outarray,"#bin";
		foreach (@bincontainer) {
			push @outarray,$_;
		}
	}
	$outarray[0] = $outarray[0] .  "\t" . "score-$group" . "\t" . "variation-$group" ;
	
	foreach $bin (@bincontainer) {
		$i++;
		$outarray[$i] = $outarray[$i] . "\t" . $binzmatrix{$bin} . "\t" . $binscoresd{$bin} ;
	}
	return (\@outarray);
	
}


sub BINZSOCRES {
	my ($chr,$bincontainer,$binmatrix,$mean_or_zscore,$group)=@_;
	my @bincontainer=@$bincontainer;
	my %binmatrix=%$binmatrix;
	my %chr = %$chr;
	my ($coor,$keys,$bin,@allbincounts,@bincounts,$allmean,$stdev,$sterr,$binmean,$binerr,$binzscore,%binzmatrix,$binsd,%binscoresd,$groupchr,$groupcoord);
	foreach ( @$group) {
		($groupchr,$groupcoord) = split(/-/,$_);
		foreach $bin ( @bincontainer ) {
			push @allbincounts,$binmatrix{$groupchr}{$groupcoord}{$bin};
		}
	}

	$allmean = MEAN (\@allbincounts); #print "mean: $allmean\n";
	($stdev,$sterr) = STDEV (\@allbincounts); #print "std: $stdev\n";
	
	foreach $bin ( @bincontainer ) {
		@bincounts=();
		foreach $keys ( keys %chr ) {
			foreach $coor (keys %$keys) {
				push @bincounts,$binmatrix{$keys}{$coor}{$bin};
				#print "$keys\t$coor\t$bin\t" , $binmatrix{$keys}{$coor}{$bin} , "\n";
			}
		}
		$binmean=MEAN (\@bincounts);
		($binsd,$binerr)=STDEV (\@bincounts);
		$binzscore=ZSCORE ( $binmean,$stdev,$allmean );
		#print "$binmean:$binzscore\t";
		if ($mean_or_zscore eq "m" ) {
			$binzmatrix{$bin}=$binmean;
			$binscoresd{$bin}=$binerr;
		} elsif ($mean_or_zscore eq "z") {
			$binzmatrix{$bin}=$binzscore;
			$binscoresd{$bin}=$binerr;
		}
	}
	#print "\n";

	return (\%binzmatrix,\%binscoresd) ;

}

sub RETURNBINCOUNTSALLCHR { 
        my($chr,$bincontainer,$exp,$control,$binwidth,$extend_length,$normalize,$totalExpReads,$totalCtrReads,$binwidth,$strand)=@_;
	my %chr = %$chr;
	my @bincontainer=@$bincontainer;
	my (%binmatrix,%chrmatrix,$keys);
	my $pm = new Parallel::ForkManager($nproc);
	$pm->run_on_finish(     # callback function to retrieve and assemble matrix from each chromosome matrices.
        sub{
                my($pid, $exit_code, $ident, $exit_signal, $core_dump, $dr) = @_;
                if(defined $dr){
			print "** $ident \tFINISHED, distribution stored. pid: $pid\n";
                        %chrmatrix = %$dr ;
                        @binmatrix{ keys %chrmatrix } = values %chrmatrix;
                }else{
                        die "Child process returned void data in 'RETURNBINCOUNTSALLCHR' with exit code: $exit_code.\n";
                }
           }
        );

	$pm->run_on_start( sub {
		my ($pid,$ident)=@_;
		print "** $ident \tstarted, calculating distribution. pid: $pid\n";
	});

	foreach $keys ( keys %chr ) {
		my $pid = $pm->start($keys) and next;
		#print "--$_\t" foreach @bincontainer ; print "\n";
		%chrmatrix = RETURNBINCOUNTS ($keys,$bincontainer,$exp,$control,$binwidth,$extend_length,$normalize,$totalExpReads,$totalCtrReads,$binwidth,$strand);
#		@binmatrix{ keys %chrmatrix } = values %chrmatrix;	
		$pm->finish(0, \%chrmatrix);
		
	}

	$pm->wait_all_children;

	return %binmatrix ;

	        foreach $bin ( @bincontainer ) {
                @bincounts=();
                foreach $keys ( keys %chr ) {
                        foreach $coor (keys %$keys) {
                                push @bincounts,$binmatrix{$keys}{$coor}{$bin};
                                print "$keys\t$coor\t$bin\t" , $binmatrix{$keys}{$coor}{$bin} , "\n";
                        }
                }}

}	

sub RETURNBINCOUNTS {
        my($chr,$bincontainer,$exp,$control,$binwidth,$extend_length,$normalize,$exptotalcounts,$controltotalcounts,$binwidth,$strand)=@_;
        #my %chr = %$chr;
	#$normalize : per cent, per thousand, per million etc. The value should be 1000000 for per million
	#$exptotalcounts : total counts in exp
	#$controltotal counts : total couns in control exp
        my @bincontainer=@$bincontainer;
		print "--$_\t" foreach @bincontainer ; print "\n";
        my %binmatrix;
        my ($coor,$keys,$bin,$expcount,$controlcount,$correctedcount,$orientation,$bincorrection);
                $keys=$chr;
                foreach $coor (keys %$keys) {
                        $orientation=$$keys{$coor} ;
                        #print "# inside $orientation";
                        foreach $bin ( @bincontainer ) {   
				#print "$bin\t";
				#if ( $bin > 0 ) {                    #
				#	$bincorrection = $bin - $binwidth ;   #
				#} else {                             #this is doneto have the start position of summit/motif to be 1 and not 0 on the plot
				#	$bincorrection = $bin ;      #the result is that correct coordinates are passed to REPORTCOUNT 
				#}                                    #
					$bincorrection = $bin ;      #the result is that correct coordinates are passed to REPORTCOUNT 
                                $expcount=REPORTCOUNTS($keys,$exp,$coor + ORIENTATION($bincorrection,$orientation) ,$binwidth,$extend_length,STRANDCORRECTION($strand,$orientation) );
				#print "-----$strand\t";
                                if ( $control ne "NULL") {
                                        $controlcount=REPORTCOUNTS($keys,$control,$coor + ORIENTATION($bincorrection,$orientation) ,$binwidth,$extend_length,STRANDCORRECTION($strand,$orientation) );
                                } else {
                                        $controlcount=0;
                                }
				if ( $normalize != 0 ) {
					$expcount = sprintf ("%.4f", $expcount * $normalize / $exptotalcounts) ;
					$controlcount = sprintf ("%.4f",$controlcount * $normalize / $exptotalcounts) ;
					#print "$expcount $controlcount\n";
				}
                                if ( $expcount > $controlcount ) {
                                        $correctedcount=$expcount-$controlcount ;
                                } else {
                                        $correctedcount=0;
                                }
                                #print "$keys\t","$coor\t$bin\t",$bin,"\t$orientation\t","$expcount\t$controlcount\t$correctedcount","\n" ;
                                $binmatrix{$keys}{$coor}{$bin}=$correctedcount;
                                #$binmatrix{$keys}{$coor}{$bin}=$correctedcount;
                        }
                        #print "\n";
                }
        return %binmatrix;
}
sub STRANDCORRECTION {
	my ($strand,$orientation)=@_;
	if ( $orientation eq "+") { return $strand ;}
	if ( $orientation eq "-" && $strand eq "fwd") { return "rev" ;}
	if ( $orientation eq "-" && $strand eq "rev") { return "fwd" ;}
	if ( $orientation eq "-" && $strand eq "both") { return "both" ;}
}


sub ORIENTATION {
        my ( $bin, $orientation )=@_;
        if ( $orientation eq "+" ) {
                return $bin ;
        } elsif ( $orientation eq "-" ) {
                return (-$bin) ;
        } else {
                return $bin ;
        }
}


sub REPORTCOUNTS {
        my ($chr,$exp,$binstart,$binwidth,$extend_length,$strand) = @_ ;
        my $coordinates ;
        my $count=0 ;
	if ($strand eq "fwd" ) { $strand = "+" } elsif ($strand eq "rev" ) { $strand = "-"; } #elsif { print "unknown strand type\n"; exit; }
        $chr=$exp.$chr;
        for $coordinates ( ( $binstart - $extend_length + 1 ) .. ( $binstart + $binwidth -1 )  ) {
                if (exists $$chr{$coordinates} ) {
			if ($strand ne "both") {
                        	$count += $$chr{$coordinates}{$strand} if exists $$chr{$coordinates}{$strand} ;
			} elsif ($strand eq "both")  {
                        	$count += $$chr{$coordinates}{"+"} if exists $$chr{$coordinates}{"+"} ;
                        	$count += $$chr{$coordinates}{"-"} if exists $$chr{$coordinates}{"-"} ;
			} else {
				print "exit: unknown strand type (should be fwd, rev or both )";
				exit;
			}
                        #print "$coordinates\texists=", $$chr{$coordinates} ,"\n" ;
                } else {
				next ;
 #                       print "$coordinates\tnope\n" ;
                }
        }

    #    print "total counts= $count \n";
        return $count ;

}

sub READSUMMITS {
        my ($summitFile,$groupSize)=@_;
        my (%chr,$orientation,@temp,$chromosome,$start,$end,$chrCoord) ;
	my $i=-1;
	my $j=0;
	my $connector="-"; 
	my $group;
	my $groupall=$summitFile."all" ;
	my @groups=($groupall);
	#push @groups,$group if $groupSize != 0 ;
        open PEAKS,$summitFile;
        while (<PEAKS>){
		chomp;
		$i++;
		if ( ( $groupSize != 0 ) && ( $i % $groupSize == 0 )  ){
			$j++;
			$group=$summitFile.$j;
			push @groups,$group;
			#print "$group $i \n";
		}
                @temp=split(/\t/,$_ ) ;
		$start=$temp[1];
		$end=$temp[2];
		$chromosome=$temp[0];
                #print "# strand $6\n";
		if ( @temp == 3 ) {
			$orientation=1;
		} elsif ( $temp[5] eq "-" || $temp[5] eq "+") {
			if ( $temp[5] eq "-"  ) {
				$start = $end - 1 ;
			}
                        $orientation=$temp[5];
                } else {
                        $orientation=1;
                }
		$chrCoord=$chromosome.$connector.$start ;
		push @$group,$chrCoord if $groupSize != 0;
		push @$groupall,$chrCoord ;
                $$chromosome{$start}=$orientation;
                $chr{$chromosome}=1;
                #print "$1\t$2\t$5\n";

        }
	close PEAKS ;
        return (\%chr,\@groups);
}


sub BINS {
        my($binwidth,$bflank,$aflank)=@_;
        my @bincontainer;
	my $i;
        for ( $i=-$bflank;$i<=$aflank;$i+=$binwidth) {
		#next if $i == 0 ;
                push @bincontainer,$i ;
        }
        return @bincontainer ;
}

sub READALIGNEDSLICEDTAGS {
	my ($slice,$bedFile,$extend_length,$rscore,$dcol) =@_ ;
	my ($totalExpReads,@temp,@split,$chr);
	@temp=@$slice;
	#open READS,$bedFile;
	#chomp (@temp=<READS>) ;
	#close READS ;
	#print STDERR "  $bedFile File Read\n";
	#while (<READS>){
	my $n=0;
	foreach (@temp) {
		next if $_ =~ /^#/ ;
        	chomp;
		@split=();
		@split=split(/\t/,$_);
		$chr=$split[0];
        	#$_=~/^(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)/;
        	$totalExpReads++;
		#print "$direction $score\n";
                $chrhash = $bedFile.$chr ; # not local variable (deliberately)
        	if ( $split[5] eq "-") {
			if ( $rscore eq "n" ) {
				$$chrhash{ $split[2] - $extend_length }{"-"}++;
			} else {
				$$chrhash{ $split[2] - $extend_length }{"-"}+=$split[$dcol - 1]  ;
			}
                	#print "$1\t",$3 - $extend_length,"\t$6\n";
        	} else {
			if ( $rscore eq "n" ){
                		$$chrhash{$split[1]}{"+"}++;
			} else {
				$$chrhash{$split[1]}{"+"}+=$split[$dcol - 1] ;
			}
                	#print "$1\t",$2,"\t$6\n";
        	}
		#print $$chrhash{$split[1]} ,"\n";
		#$chr{$1}=1;
       		#print "$chrhash\t$split[0] $split[1] $split[2] $split[3] $split[4] $split[5] \n";
		$n++;
		if ($n % 1000000 == 0 ){
			print STDERR "  processing $bedFile : reads processed $n\n";
		}
	}
	print STDERR "  extending done\n";
	return $totalExpReads ;
}

sub READALIGNEDTAGS {
        my ($bedFile,$extend_length,$rscore,$dcol) =@_ ;
        my ($totalExpReads,@temp,@split,$chr,$each);
        #@temp=@$slice;
        open READS,$bedFile;
        chomp (@temp=<READS>) ;
        close READS ;
        print STDERR "  $bedFile File Read\n";
        #while (<READS>){
        my $n=0;
        #foreach (@temp) {
	while ( $each = shift @temp ) {
                next if $each =~ /^#/ ;
                chomp;
                @split=();
                @split=split(/\t/,$each);
                $chr=$split[0];
                #$_=~/^(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)/;
                $totalExpReads++;
                #print "$direction $score\n";
                $chrhash = $bedFile.$chr ; # not local variable (deliberately)
                if ( $split[5] eq "-") {
                        if ( $rscore eq "n" ) {
                                $$chrhash{ $split[2] - $extend_length }{"-"}++;
                        } else {
                                $$chrhash{ $split[2] - $extend_length }{"-"}+=$split[$dcol - 1]  ;
                        }
                        #print "$1\t",$3 - $extend_length,"\t$6\n";
                } else {
                        if ( $rscore eq "n" ){
                                $$chrhash{$split[1]}{"+"}++;
                        } else {
                                $$chrhash{$split[1]}{"+"}+=$split[$dcol - 1] ;
                        }
                        #print "$1\t",$2,"\t$6\n";
                }
                #print $$chrhash{$split[1]} ,"\n";
                #$chr{$1}=1;
                #print "$chrhash\t$split[0] $split[1] $split[2] $split[3] $split[4] $split[5] \n";
                $n++;
                if ($n % 1000000 == 0 ){
                        print STDERR "  processing $bedFile : reads processed $n\n";
                }
        }
        print STDERR "  extending done\n";
        return $totalExpReads ;
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

sub PRINTTIME {
	my $msg = shift @_;
	my $time = time - $time_tag ;
	my $total_time = time - $start_time ;
	$time = convert_time($time);
	$total_time = convert_time($total_time);
	print ("Total time elapsed:$total_time\t $time taken in $msg\n\n");
	$time_tag = time;
}

sub convert_time {
  my $time = shift;
  my $days = int($time / 86400);
  $time -= ($days * 86400);
  my $hours = int($time / 3600);
  $time -= ($hours * 3600);
  my $minutes = int($time / 60);
  my $seconds = $time % 60;

  $days = $days < 1 ? '' : $days .'d ';
  $hours = $hours < 1 ? '' : $hours .'h ';
  $minutes = $minutes < 1 ? '' : $minutes . 'm ';
  $time = $days . $hours . $minutes . $seconds . 's';
  return $time;
}
