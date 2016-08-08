#!/usr/bin/perl 

($feature_file)=@ARGV;

open OPEN,$feature_file;
chomp(@features=<OPEN>);

$order_column=15;

for ($i=5;$i<=50;$i+=5){
	&SVM(\@features, 15 , $i );
}

sub SVM {
	my($file,$order_column,$pc)=@_;
	my @fileoriginalorder=@$file;
	my $total=@fileoriginalorder;
	my $head_stop=sprintf ("%.0f", $pc * $total / 100) ;
	my $tail_stop=$total - $head_stop; 
	@file = sort { (split /\t/,$b)[$order_column] <=> (split /\t/,$a)[$order_column] } @fileoriginalorder ;
	#print "$order_column \t $pc\t $total\t $head_stop \t$tail_stop \n";
	my @positive=@file[0 .. $head_stop] ;
	my @negative=@file[$tail_stop .. $#file];
	foreach ( 0 .. $#negative) {
		#print "$negative[$_]\n";
		$negative[$_] =~ s/\t/ /g;
		$negative[$_] =~ s/^1/-1/;
	}
	foreach (0 .. $#positive) {
		$positive[$_] =~ s/\t/ /g;
	}
	open POSWRITE, ">$feature_file.$pc.features.txt";
	print POSWRITE "$_\n" foreach @positive ;

	#open NEGWRITE, ">$feature_file.$pc.neg.txt";
	print POSWRITE "$_\n" foreach @negative ;
	
	system ( "svm_learn $feature_file.$pc.features.txt $feature_file.$pc.model.txt");
	system ("svm_classify $feature_file $feature_file.$pc.model.txt $feature_file.$pc.model.pred");
	open PRED, "$feature_file.$pc.model.pred";
	chomp (my @predicted=<PRED>);
	open WRITEPRED,">$feature_file.$pc.model.pred.out";
	foreach ( 0 .. $#fileoriginalorder){
		if ( $predicted[$_] <= 0 ) {next;} else {
			$fileoriginalorder[$_]=~s/^-1|^1/$predicted[$_]/;
			print WRITEPRED $fileoriginalorder[$_] , "\n" ,
		}
	}
	print "Percentage: $pc\n";
	&REPORT_COVERAGE("$feature_file.$pc.model.pred.out", 27341 ) ;
	print "\n\n-next cycle\n";
}



#REPORT_COVERAGE("see.txt",10551);

sub REPORT_COVERAGE {
	my ($predictions,$total)=@_;
	my %covHash;
	open OPEN,$predictions;

	while (<OPEN>){
		$covHash{(split /\t/,$_)[6]}++;
	}
	
	#print "$_\t$covHash{$_}\n" foreach keys %covHash ;
	my @covered=();
	my $sumOfPredictedMotifs=0;
	foreach (keys %covHash) {
		#print $_ , "\t" ,$covHash{$_} , "\n";
		$sumOfPredictedMotifs += $covHash{$_} ;
		push @covered,$covHash{$_};

	}
	my $mean=MEAN(\@covered);
	my $total_predicted=$#covered + 1 ;
	my $sd=STDEV(\@covered); 
	my $coverage=sprintf ("%.2f",@covered * 100 / $total);
	print "totalPredicted=$total_predicted\t$sumOfPredictedMotifs\tmean=$mean\tsd=$sd\tcoverage=$coverage\n";
	
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

