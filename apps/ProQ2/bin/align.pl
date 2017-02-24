#!/usr/bin/perl -w
#use lib '/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/bioperl-live/';
#use lib '/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/bioperl-live/bioperl-ext-1.5.1/';
use lib '/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/bioperl-live/bioperl-ext-1.5.1/Bio/Ext/blib/lib/';
use lib '/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/bioperl-live/bioperl-ext-1.5.1/Bio/Ext/blib/arch/';
#BEGIN {
## #   push(@INC,'/afs/pdc.kth.se/home/a/arnee/MODULES/perl5/lib/site_perl/5.6.0/i386-linux/');
##   # push(@INC,'/opt/local/lib/perl5/site_perl/5.8.8/darwin-2level/auto/');
#    unshift(@INC,'/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/bioperl-live/bioperl-ext-1.5.1/');
#}

use Bio::Ext::Align;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Tools::pSW;
use Bio::LocatableSeq;
use Bio::Seq;

#BEGIN {
# #   push(@INC,'/afs/pdc.kth.se/home/a/arnee/MODULES/perl5/lib/site_perl/5.6.0/i386-linux/');
#   # push(@INC,'/opt/local/lib/perl5/site_perl/5.8.8/darwin-2level/auto/');
#    unshift(@INC,'/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/bioperl-live/');
#}

sub align   # Takes two strings removes all dashes and returns the alignment. 
{
    my ($seq1,$seq2)=@_;
    my $input1=$seq1;
    my $input2=$seq2;
    $seq1=~s/-//g;
    $seq2=~s/-//g;
    #$seq1=remove_dashes($seq1);
    #$seq2=remove_dashes($seq2);
    $seq1=~s/\n//g;
    $seq2=~s/\n//g;
    $seq1=~s/\s+//g;
    $seq2=~s/\s+//g;
   # if(length($seq1)==1 || length($seq2)==1)
   # {
   #	my $len1=length($seq1); 
   #	my $len2=length($seq2);
   #	my $len3=length($input1);
   #	my $len4=length($input2);
   #	if($len1==1)
   #	{
   #	   $dashes=length($input2)-length($input1);
   #	   $ali_return1=$input1._dashes($dashes);
   #	   $ali_return2=$input2;
   #	}
   #	if($len2==1)
   #	{
   #	   $dashes=length($input1)-length($input2);
   #	   $ali_return1=$input1;
   #	   $ali_return2=$input2._dashes($dashes);
   #	}
   #	return ($ali_return1,$ali_return2)
   # }
   #
   #
   # else
   # {
    my ($ali_return1,$ali_return2);
    my $factory=new Bio::Tools::pSW('-matrix' => '/afs/pdc.kth.se/home/b/bjornw/bjorn/bioperl/bioperl-live/examples/blosum62.bla','-gap' => 1,'-ext' => 0);
    my $seq_obj1=Bio::Seq->new(-moltype => 'protein', -seq => $seq1, -id => "seq1");
    my $seq_obj2=Bio::Seq->new(-moltype => 'protein', -seq => $seq2, -id => "seq2");
    my $aln = $factory->pairwise_alignment($seq_obj1,$seq_obj2);
    #my $alnout = new Bio::AlignIO(-format => 'fasta',
     #                             -fh     => \*STDOUT);

    #$alnout->write_aln($aln);
#    print "align\n";
#    print $aln."\n";
#    foreach my $tmp(keys(%{$aln}))
#    {
#	print $tmp."\n";
#	print $aln{$tmp};
   # }
   # print $seq_obj2->seq();
   # print "\n";
   # print $aln,"\n";
    #my $nice_ali=$factory->align_and_show($seq_obj1,$seq_obj2,*STDOUT);
    
    ($ali_return1,$ali_return2)=fix_alignment($aln,$seq_obj1,$seq_obj2);
    #print "$ali_return1\n$ali_return2\n";
    return ($ali_return1,$ali_return2);
    
}

sub fix_alignment    # $seq1 and $seq2 must be in the same order as in the $aln otherwise it gets wrong.
{
    my($aln,$seq1,$seq2)=@_;
    my $ali_seq1="";
    my $ali_seq2="";
    # Parse the alignment
    my ($seq,@start, @end,@ali_seqs);
    my $i=0;

    foreach $seq ($aln->each_seq())
    {
	#print $seq->seq()."\n";
	$start[$i]=$seq->start();
	$end[$i]=$seq->end();
	$ali_seqs[$i]=$seq->seq();
	#print $seq."\n";
	#print $i."\n";
	#print $ali_seqs[$i]."\n";

	#print "$start[$i] $end[$i]\n";
	$i++;
    }
    #print $start[0]."\n".$start[1]."\n";
# Reformat alignment so that it contain all resides and -.
    
  #  print "$ali_seqs[0]\n$ali_seqs[1]\n";
  #  exit;
# Fix the begining
    if($start[0]!=1)
    {
	$ali_seq1.=$seq1->subseq(1,$start[0]-1);
	$ali_seq2.=_dashes($start[0]-1);
    }
    if($start[1]!=1)
    {
	$ali_seq1.=_dashes($start[1]-1);
	$ali_seq2.=$seq2->subseq(1,$start[1]-1);
    }
    
# Add the alignment
    
    $ali_seq1.=$ali_seqs[0];
    $ali_seq2.=$ali_seqs[1];
    
# Fix alignment end;
 #   print "$end[0] ".$seq1->length()."\n";
    if($end[0]<$seq1->length())
    {
	my $len=$seq1->length();
	$ali_seq1.=$seq1->subseq($end[0]+1,$len);
	$ali_seq2.=_dashes($len-$end[0]);
    }
    #print "$end[1] ".$seq2->length()."\n";
    if($end[1]<$seq2->length())
    {
	my $len=$seq2->length();
	$ali_seq1.=_dashes($len-$end[1]);
	$ali_seq2.=$seq2->subseq($end[1]+1,$len);
    }
    return ($ali_seq1,$ali_seq2);
}

1;


sub _dashes
{
    my $number=shift;
    my $str="";

    for(my $i=0;$i<$number;$i++)
    {
	$str.="-";
    }
    return $str;

}

sub merge_ali   # $A aligned with $B1, $B2 aligned with $C, returns $A aligned with $C.
{               # $B1 and $B2 have exactly the same residues.
    my($A,$B1,$B2,$C) = @_;
    my $A_return="";
    my $C_return="";
    my @A=split('',$A);
    my @B1=split('',$B1);
    my @B2=split('',$B2);
    my @C=split('',$C);
    my $len1 = scalar @A; #length of first aligment pair 
    my $len2 = scalar @B2; #length of second aligment pair

    my $j=0;   # $i indexes in the first alignment $j in the second
    my $i=0;
    my $temp="";;
    my $temp2="";
    while($i<$len1 || $j<$len2) 
    {
	#print "$B1[$i] $B2[$j]\n";
	if(defined($B1[$i]) && $B1[$i] eq "-")
	{
	    $temp2.="-";
	    $temp.=$A[$i];
	    $i++;
	}
	elsif($B2[$j] eq "-")
	{
	    $temp2.=$C[$j];
	    $temp.="-";
	    $j++;
	}
	elsif($B1[$i] eq $B2[$j])
	{
	    #print "$i $j $A[$i] $C[$j]\n";
	    $temp2.=$C[$j];
	    $temp.=$A[$i];
	    $i++;
	    $j++;
	}
	
    }
    my $temp3="";
    my $temp4="";
    my @temp=split('',$temp);
    my @temp2=split('',$temp2);
    #print $temp."\n\n".$temp2."\n\n";
    # Remove all positions which has gaps in both sequences.
    for(my $i=0;$i < scalar @temp2;$i++)
    {
	if(not ($temp[$i] eq "-" && $temp2[$i] eq "-"))
	{
	    $temp3.=$temp[$i];
	    $temp4.=$temp2[$i];
	}

    }

    #print $temp3."\n\n".$temp4."\n\n";
    return($temp3,$temp4);
}
