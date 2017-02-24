#!/usr/bin/perl -w


my $caspdir="/afs/pdc.kth.se/home/b/bjornw/CASP9/";
my $author="1239-5930-4200";

foreach my $seqfile(glob("$caspdir/T0*/sequence"))
{
    if($seqfile=~/(T0\d\d\d)/) {
	my $target=$1;

	my $outdir=dirname($seqfile);
	my @list=`grep Global $outdir/models/*.proq2|sort -n -k -r 3`;
	foreach my $line(@list) {

	    print "$target $line\n";
	}



    }

}


