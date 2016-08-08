#!/usr/bin/perl 

$0 =~ /(.*)\/(.*)/;
$prog_dir=$1;

($feature_svm_outfile,$negative_features)=@ARGV;

#./bin/eval_roc.pl VEZF1_macs_subpeaks_VEZF1_6848_motif_10_20/motif_master_features.bed.15.model.pred.out null_VEZF1_macs_subpeaks_VEZF1_6848_motif_10_20/motif_master_features.bed VEZF1_macs_subpeaks_VEZF1_6848_motif_10_20/motif_master_features.bed.15.model.txt

$feature_svm_outfile=~/(.*)\/(.*)/;
$dir=$1;
$file=$2;

print "file\tpercent\tAccuracy\tprecision\trecall\n";
for ($i=5;$i<=50;$i+=5){
        #&SVM(\@features, 15 , $i );
	$predicted_positives_features="$dir/$file.$i.model.pred.out";
	$model_file="$dir/$file.$i.model.txt";
	#print "$predicted_positives_features\t$negative_features\t$model_file\n";
	my ($accuracy,$correct,$incorrect,$total,$precision,$recall)=MAKErun($predicted_positives_features,$negative_features,$model_file);
	system ("Rscript $prog_dir/plotROC.R $dir/roc/$file.$i.model.pred.out.class $dir/roc/$file.$i.model.pred.out.class.out  > $dir/roc/$file.$i.model.pred.out.class.out.rout");
	print "$predicted_positives_features\t$i\t$accuracy\t$precision\t$recall\n";

}
#exit;

($accuracy,$correct,$incorrect,$total,$precision,$recall)=MAKErun($feature_svm_outfile,$negative_features,$modelfile);
print "acc:$accuracy\tprec:$precision\trecall:$recall\n";

sub MAKErun {
	my ($feature_svm_outfile,$negative_features,$modelfile)=@_;
	$feature_svm_outfile=~/(.*)\/(.*)/;
	my $dir=$1;
	my $file=$2;
	#print "$dir $file";

	if (! -d "$dir/roc" ) {
		mkdir "$dir/roc";
	}

	my @positives=getPositives($feature_svm_outfile);
	#$totalpositives=$#positives;
	#print "total pos = $totalpositives\n";
	my @negatives=getNegatives($negative_features,$#positives,15);

	open WRITE, ">$dir/roc/$file.class";
	print WRITE "$_\n" foreach @positives;
	print WRITE "$_\n" foreach @negatives;

	system ("svm_classify $dir/roc/$file.class $modelfile $dir/roc/$file.class.out $modelfile > $dir/roc/$file.class.log");
	open OPEN,"$dir/roc/$file.class.log";
	my ($accuracy,$correct,$incorrect,$total,$precision,$recall);
	while (<OPEN>){
		if (/^Accuracy on test set: (.*?)% \((\d+) correct, (\d+) incorrect, (\d+) total\)/){
			$accuracy=$1;
			$correct=$2;
			$incorrect=$3;
			$total=$4;
		}	
		if (/^Precision\/recall on test set: (.*?)%\/(.*?)%/){
			$precision=$1;
			$recall=$2;
		}
	}

#Accuracy on test set: 97.20% (261554 correct, 7530 incorrect, 269084 total)
#Precision/recall on test set: 94.97%/100.00%
	return ($accuracy,$correct,$incorrect,$total,$precision,$recall);
}


sub getNegatives {
	my($feature_svm_outfile,$number,$order_column)=@_;
	my (@negatives,$i);
	open OPEN,$feature_svm_outfile;
	chomp(my @feature_svm_outfile=<OPEN>);
	@feature_svm_outfile_sorted_by_cov = sort { (split /\t/,$b)[$order_column] <=> (split /\t/,$a)[$order_column] } @feature_svm_outfile ;
	foreach (@feature_svm_outfile_sorted_by_cov) {
		$i++;
		my @split=split(/\t/,$_);
		$split[0]= -1;
		my $join=join "\t",@split;
		push @negatives,$join;
		last if ( $i == $number);
	}
	return (@negatives)
}


sub getPositives {
	my($feature_svm_outfile)=@_;
	my @positives;
	open OPEN,$feature_svm_outfile;
	chomp (my @feature_svm_outfile=<OPEN>);
	foreach (@feature_svm_outfile){
		my @split=split(/\t/,$_);
		$split[0]=1;
		#print $split[0] , "\n";
		my $join=join "\t",@split;
		push @positives,$join;
	}
	return (@positives);
}
