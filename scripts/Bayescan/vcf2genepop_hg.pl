#!/usr/bin/perl

$usage= "

vcf2genepop.pl :

Converts VCF to multiallelic GENEPOP, preserves chromosome and position info

Arguments: 
vcf=[filename] : vcf file to convert
   pops=[list] : comma-separated perl patterns to determine population 
                 affiliation of samples in the VCF file based on sample names

Output:
GENEPOP formatted genotypes, printed to STDOUT

Example:
vcf2genepop.pl vcf=filtered.vcf pops=O,K,M,S > filtered.gen

";

%dict = ('A', 1, 'C', 2, 'G', 3, 'T', 4);

my $vcf="";
my @pops=();

if ("@ARGV"=~/vcf=(\S+)/) {$vcf=$1;}
if ("@ARGV"=~/pops=(\S+)/) { @pops=split(/\,/, $1);}

open VCF, $vcf or die "cannot open input $vcf\n\n$usage";

my @samples=();
my $header;
while (<VCF>) {
	chomp;
	if ($_=~/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t(.+)/) {
		@samples=split(/\s/,$1);
		last;
	}
}
my $nsam=$#samples+1;
my @loci=();
my %lgts={};
my @gts=();
my $idx=0;
while (<VCF>) {
	chomp;
	my $line=$_;
	if ($line=~/(\S+)\t(\S+)\t\S+\t(\S+)\t(\S+)\t\S+\t\S+\t\S+\t\S+\t(.+)/ ){
		$idx = $idx + 1;
		$REF = $3;
		$ATL = $4;
		@swaps = ($REF,$ATL);

	
		
		my $locus=$1."_".$2;
		
		push @loci, $locus;
		@gts=split(/\s/,$5);
		for($g=0; $gt=$gts[$g];$g++) {
			($gtt,@rest)=split(/:/,$gt);
			($gt1,$gt2)=split(/[|\/]/,$gtt);

			if ($gt1 eq ".") {
				$gt1=0;
				$gt2=0;
			}
			else {
				#print @swaps,'---',$gt1;
				$gt1 = $swaps[$gt1];
				#print '-->',$gt1,"\n";
				$gt1 = $dict{$gt1};
				
				$gt2 = $swaps[$gt2];
				#print %dict,'---',$gt2;
				$gt2 = $dict{$gt2};
				#print '-->',$gt2,"\n";

				#@previous design
				#$gt1=$gt1+1;
				#$gt2=$gt2+1;
			}

			$lgts{$locus}{$samples[$g]}="0".$gt1."0".$gt2;
			print $lgts
		}
	}
}

@samples=sort(@samples);

print "converted from $vcf\n";
print join(",\n",@loci), "\n";
my $newpop;
my $oldpop;

foreach $sa (@samples){
	foreach $p (@pops) {
	        if ($sa=~/$p/) {$newpop=$p;}
	}
	if ($newpop ne $oldpop) {
	        print "POP\n";
	       	$oldpop=$newpop;
	}
	print $sa, ",";
	foreach $l (@loci) {
		print " ", $lgts{$l}{$sa};
	}
	print "\n";
}
