#!/usr/bin/perl 
($file,$cutoff)=@ARGV;

($data,$initialpeaks)=parseAllData($file);
#print "$totalpeaks\n";
#print "$_\n" foreach @$data;

if (! -d "$dir/roc" ) {
	mkdir "$file.distribution";
} else {
	print "Directory $file.distribution already exists\n";
}


#$log=negativelog10(5e-05);

#print "log-- $log\n";
#$filter=getFiltered($data,$cutoff);
#print "$_\n" foreach @$filter;
#($motifs,$mean,$STDEV)=countRANKS($filter);
#print "mean = $mean +- $STDEV\n";
#print "$_\n" foreach @$motifs ;

for ($i=0;$i<=7;$i+=0.5){
	#print "$i\n";
	my $filter=getFiltered($data,$i);
	if (@$filter==0){
		print "$i\t0\t0\t0\n";
		next;
	}
	my ($motifs,$mean,$STDEV,$totalpeaks)=countRANKS($filter);
	$totalpeakspc= sprintf("%.1f", $totalpeaks * 100 / $initialpeaks );
	open WRITE,">$file.distribution/$file.distribution.$i.motifs.txt";
	print  WRITE "$_\n" foreach @$motifs ;
	print "$i\t$totalpeakspc\t$mean\t$STDEV\t$totalpeaks\n";
	system ("weblogo -f $file.distribution/$file.distribution.$i.motifs.txt -o $file.distribution/$file.distribution.$i.motifs.png -F png");
}


sub countRANKS {
	my ($filtered)=@_;
	my @filtered=@$filtered;
	my (@distances,@motifs,%ranks);
	foreach (@filtered){
		my @split=split(/\t/,$_);
		$ranks{ $split[0] } = 1;
		push @distances,$split[3] ;
		push @motifs,$split[2];
	}
	my $totalpeaks=scalar(keys %ranks);
	my $mean=sprintf("%.1f",MEAN(\@distances));
	my $STDEV=sprintf("%.1f",STDEV(\@distances));
	return (\@motifs,$mean,$STDEV,$totalpeaks);
}



sub getFiltered {
	my ($all,$cutoff)=@_;
	my @filter;
	my @all=@$all;
	foreach (@all){
		my @split=split(/\t/,$_);
		if ( ( negativelog10 ( $split[1] ) ) > $cutoff ) {
			push @filter,$_;
		}
	}
	return (\@filter);
}



sub parseAllData {
	my ($file)=@_;
	my @data;
	my %allranks;
	my $totalpeaks;
	open OPEN,$file;
	while (<OPEN>){
		next if $_=~/^#/ ;
		my @split=split(/\t/,$_);
		my $rank=$split[3];
		my $evalue=$split[4];
		my $motif=$split[6];
		my $summit=$split[12];
		my $start=$split[1];
		my $end=$split[2];
		my $motifdistfromsummit=  abs int ( $start + (($end - $start)/2) - $summit ) ;
		$allranks{$rank}=1;
		push @data, "$rank\t$evalue\t$motif\t$motifdistfromsummit";
	}
	$totalpeaks=scalar (keys %allranks);
	return (\@data,$totalpeaks);
}

sub negativelog10 {
	my $value = shift;
	$log10 = log10($value);
	return (-$log10);
}

sub log10 {
	my $n = shift;
	return log($n)/log(10);
}

sub MEAN{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
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
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

