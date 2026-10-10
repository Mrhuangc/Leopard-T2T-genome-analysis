use warnings;
use strict;

my $gtf = shift or die "Usage:\n    perl $0 gft > gff3\n";

open(F,$gtf) || die "Can't open the $gtf file!\n";
my $gene;
my $mrna;
my $cds_save;
my %count;
print "##gff-version 3\n";
while(<F>){
	chomp;
	next if(/^#|^\s*$/);
	my @a = split /\t/;
	if($a[2] eq "gene"){
		$gene = $a[8];
		$a[8] = "ID=$gene;Name=$gene";
		print "\n" if($. != 1);
		print join("\t",@a[0..8]),"\n";
	}elsif($a[2] eq "transcript"){
		$a[2] = "mRNA";
		$mrna = $a[8];
		$a[8] = "ID=$mrna;Parent=$gene";
		print join("\t",@a[0..8]),"\n";
	}elsif($a[2] eq "CDS"){
		$count{$gene}{$mrna}{cds}++;
		$a[8] = "ID=CDS-${mrna}.$count{$gene}{$mrna}{cds};Parent=$mrna";
		$cds_save = join("\t",@a[0..8])."\n";
	}elsif($a[2] eq "exon"){
		$count{$gene}{$mrna}{exon}++;
		$a[8] = "ID=exon-${mrna}.$count{$gene}{$mrna}{exon};Parent=$mrna";
		print join("\t",@a[0..8]),"\n";
		print $cds_save if($cds_save);
	}else{
		next;
	}
}
close F;
