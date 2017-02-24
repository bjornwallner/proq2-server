#!/usr/bin/perl -w 
require '/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/ProQ2/bin/align.pl';
$subunits=1;
$ss2=$ARGV[0];
$seqfile=$ARGV[1];
$outfile=$ARGV[2];
$subunits=$ARGV[3] if(defined($ARGV[3]));
$seq=`grep -v '>' $seqfile`;
$seq=~s/\n//;

$profile_seq=`cat $ss2|head -n 2 |tail -n 1`;
$ss_pred=`cat $ss2 |tail -n 1`;


#print $seq."\n";
#print $profile_seq."\n";
#print $ss_pred."\n";
chomp($profile_seq);
chomp($ss_pred);
@rsa=split(//,$ss_pred);

my $temp_seq="";
my @temp_rsa=();
for(my $i=0;$i<$subunits;$i++) {
    $temp_seq.=$profile_seq;
    @temp_rsa=(@temp_rsa,@rsa);
}
$profile_seq=$temp_seq;
@rsa=@temp_rsa;


($a1,$a2)=align($profile_seq,$seq);

print "P: $a1\nS: $a2\n";
#exit;
@a1=split(//,$a1);
@a2=split(//,$a2);

#$outfile=$seqfile;
#$outfile=~s/\.seq/\.acc/g;
#print $outfile."\n";
open(OUT,">$outfile");
#print OUT "$seq\n";
$resno=0;

foreach $rsa(@rsa) {
    if($a1[$resno] eq $a2[$resno]) {
	    
	print OUT $rsa;
    }
    $resno++;


}
print OUT "\n";


