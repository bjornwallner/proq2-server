#!/usr/bin/perl -w 
require '/local/www/services/ProQ2/apps/ProQ2/bin/bjornlib.pl';
#OBS variable names should be changed 
$subunits=1;
$ss2=$ARGV[0];
#$seqfile_acc=$ARGV[1];
$seqfile=$ARGV[1];
$outfile=$ARGV[2];

$subunits=$ARGV[3] if(defined($ARGV[3]));
$seq=`grep -v '>' $seqfile`;
$seq=~s/\n//;

#$profile_seq=`grep -v '>' $seqfile_acc`;
#$profile_seq=~s/\n//;

$profile_seq=`cat $ss2|head -n 2 |tail -n 1`;

$ss_pred=`cat $ss2 |tail -n 1`;
$ss_pred=~s/\n//;


#print "$seqfile: $seq\n";
#print "$seqfile_model: $profile_seq\n";
#print $ss_pred."\n";
#print "=====\n";
#exit;
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


