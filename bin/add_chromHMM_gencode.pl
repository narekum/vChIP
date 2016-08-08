#!/usr/bin/perl 
($state,$master)=@ARGV;

open OPEN,$state;
while (<OPEN>){
	chomp;
	@split=split (/\t/,$_);
	$hmmlocation="$split[6]:$split[7]-$split[8]" ;
	$tss_location="$split[16]:$split[17]-$split[18]";
	#$gene_id "ENSG00000227232.4"; gene_type "pseudogene";
	$split[22] =~ /gene_id "(.*?); gene_type "(.*?)";/ ; 
	$gene_id=$1;
	$gene_type=$2;
	$p=join "\t",$hmmlocation,$split[9],$tss_location,$gene_id,$gene_type,$split[23];
	#print "$p\n" ;
	$hmmstate_tss{ $split[3] } = $p ;
}
#chr1:28337-29737	1_Active_Promoter	chr1:29369-29370	ENSG00000227232.4"	pseudogene	14
open OPEN,$master;
while (<OPEN>){
	chomp;
	if ($_=~/^#/) {
		print "$_\tChrom_HMM_Location\tChrom_HMM_state\tGencode_TSS_location\tGencode_gene_id\tGencode_gene_type\tDistance_from_TSS\n";
	} else {
		@split=split (/\t/,$_);
		print "$_\t$hmmstate_tss{ $split[3] }\n";
	}
	
}
