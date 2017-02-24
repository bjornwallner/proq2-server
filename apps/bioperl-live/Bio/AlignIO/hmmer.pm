#
# BioPerl module for Bio::AlignIO::casp
#
#	
#
# You may distribute this module under the same terms as perl itself
# _history
# September 5, 2000
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::hmmer - hmmer sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::AlignIO class.

=head1 DESCRIPTION

This object can transform Bio::SimpleAlign objects to and from fasta flat
file databases.

=head1 FEEDBACK

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHORS - Björn Wallner

Email: bjorn@sbc.su.se


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


package Bio::AlignIO::hmmer;
use vars qw(@ISA);
use strict;

use Bio::AlignIO;

@ISA = qw(Bio::AlignIO);


=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
 Returns : SimpleAlign object - returns 0 on end of file
	    or on error
 Args    : NONE

=cut


#d/d1jw3a_.d.208.1.1: domain 1 of 1, from 5 to 139: score 233.5, E = 3.6e-67
sub next_aln {
    my $self = shift;
    my $line;
    my $parent;
    my @target_ali=();
    my @template_ali=();
    my $target_ali="";
    my $template_ali="";
    
    my $aln =  Bio::SimpleAlign->new();
    my $jump_to_next;
    my $header="";
    my $first_aligned=0;
    my $last_aligned=0;
    my $method;
    my ($start,$end,$score,$e_val,$template_id);
    my $name="undef";
    my $get_target=0;
    while(defined ($line = $self->_readline)) # This will read in the one line of file.
    {
	chomp($line);
	#print $line,"\n";
	if($line=~/\w/ && $get_target==1)
	{
	    #print "$line\n";
	    $line=~s/\*\->//g;
	    $line=~s/<\-\*//g;
	    $line=~s/\s+//g;
	    
	    $target_ali.=uc($line);
	    #print "TARGET $target_ali\n";
	    
	    $get_target=0;
	}

	if($line=~/$name/) #template_ali
	{
	    #print "$line\n";
	    my $temp=substr($line,13,62);
	    $temp=~s/\d+//g;
	    $temp=~s/\s+//g;
	    #print "TEMPLATE $temp\n";
	    $template_ali.=uc($temp);
	    #print "TEMPLATE $template_ali\n";
	    if($line=~/$end\s+$/) #RETURN ALIGNMENT
	    {
		$get_target=0;
		#print "LAST $line\n";
		$target_ali=~s/\./\-/g;
	#	print "\n\n",$target_ali,"\n",$template_ali,"\n";
		my $header="REMARK PARENT $template_id $parent¤";
		$header.="REMARK SCORE $e_val $score¤";
		$header.="REMARK METHOD hmmer¤";
		my $target_seq = new Bio::LocatableSeq(-moltype => "protein",
						   -seq => $target_ali,
						   -id => "Target: ");
		my $template_seq = new Bio::LocatableSeq(-moltype => "protein",
						     -seq=> $template_ali,
						     -id=> "$header");
		$aln->add_seq($target_seq);
		$aln->add_seq($template_seq); 
		return $aln;
	    }
	    $get_target=1;
	}
	if($line=~/(.+):\s+domain\s\d+\sof\s\d+,\sfrom\s(\d+)\sto\s(\d+):\sscore\s([\d\.\-]+),\sE\s=\s([\d\.e\-\+]+)/)
	#if($line=~/(.+):\sdomain\s\d+\sof\s\d+,/)
	{
	    #print "HERE: $line\n";
	    #my @temp=split(/\./,$1);
	    $parent=$1;
	    $name=substr($parent,0,10);
	    $template_id=uc(substr($parent,3,6));
	    $start=$2;
	    $end=$3;
	    $score=$4;
	    $e_val=$5;
	    $get_target=1;
	    #print "$name $start $end $score $e_val\n";
	  #  print "E-value $5\n";
	}

	    

	#if($line=~/^REMARK PARENT/)
	#{
	#    $jump_to_next=0;
	#    #my @temp=split(/\s+/,$line);
	 #   #$parent=$temp[2];
	 #   #print "PARENT:::::::::::::::::: $parent\n";
	 #   $header.=$line."¤";
	#}
	#if($line=~/^REMARK LATE ERROR/ || ($line=~/^TER/ && scalar(@target_ali) == 0))
	#{
	#    $jump_to_next=1;
	#}
	#elsif($line=~/^REMARK ALIGNMENT METHOD/)
	#{
	#    $jump_to_next=0;
	#    #$header.=$line."¤";
	#    my @temp=split(/REMARK ALIGNMENT METHOD:\s+/,$line);
	#    $method=$temp[1];
	#    $method="pp" if($method eq "PROFILE PROFILE ALIGNMENTS");
	#    $method="pp2" if($method eq "PROFILE PROFILE ALIGNMENTS 2");
	#    $method="ps" if($method eq "PROFILE SEQUENCE ALIGNMENTS");
	#    $method="ps2" if($method eq "PROFILE SEQUENCE ALIGNMENTS 2");
	#    $method="psss" if($method eq "PROFILE SEQUENCE ALIGNMENTS WITH SS");
	#    $method="psss2" if($method eq "PROFILE SEQUENCE ALIGNMENTS WITH SS 2");
	#    $method="sp" if($method eq "SEQUENCE PROFILE ALIGNMENTS");
	#    $method="sp2" if($method eq "SEQUENCE PROFILE ALIGNMENTS 2");
	#    $method="spss" if($method eq "SEQUENCE PROFILE ALIGNMENTS WITH SS");
	#    $method="spss2" if($method eq "SEQUENCE PROFILE ALIGNMENTS WITH SS 2");
	#    $method="gpp" if($method eq "GLOBAL PROFILE PROFILE ALIGNMENTS");
	#    $method="gps" if($method eq "GLOBAL PROFILE SEQUENCE ALIGNMENTS");
	#    $method="gpsss" if($method eq "GLOBAL PROFILE SEQUENCE ALIGNMENTS WITH SS");
	#    $method="gpsss2" if($method eq "GLOBAL PROFILE SEQUENCE ALIGNMENTS WITH SS 2");
	#    $method="gsp" if($method eq "GLOBAL SEQUENCE PROFILE ALIGNMENTS");
	#    $method="gspss" if($method eq "GLOBAL SEQUENCE PROFILE ALIGNMENTS WITH SS");
	#    $method="gspss2" if($method eq "GLOBAL SEQUENCE PROFILE ALIGNMENTS WITH SS 2");
	#    $header.="REMARK ALIGNMENT METHOD: $method"."¤";
	#}
	#elsif($line=~/^REMARK/ || $line=~/^SCORE/)
	#{
	#    $header.=$line."¤";
	#}
	#elsif($line=~/^\s+(\w)\s+(\d+)\s+(\w)\s+(\d+)/ && not $jump_to_next) # Does the line begin with one character?
	#{
#	#    print "test $1 $2 $3 $4\n";
	#    if ($2 > 0 && $4>0){
	#	if($method eq "ps" || $method eq "psss" || $method eq "ps2" || $method eq "psss2" || $method eq "gps" || $method eq "gpsss" || $method eq "gpsss2")  # switch target and template alignments
	#	{
	#	    my $index=$4-1;
	#	    $target_ali[$index]=$3;
	#	    $template_ali[$index]=$1;
	#	}
	#	else
	#	{
	#	    my $index=$2-1;
	#	    $target_ali[$index]=$1;
	#	    $template_ali[$index]=$3;
	#	    
	#	}
	#    #All . indicating gaps are ignored.
	#    #my $index=$2-1;
	#    #$target_ali[$index]=$1;
	#    #$template_ali[$index]=$3;
	#    }
	#}
	#elsif($line=~/^TER/ && not($jump_to_next))
	#{
	#    for(my $i=0;$i<scalar @target_ali;$i++)
	#    {
	#	if(not(defined($target_ali[$i])))
	#	{
	#	    $target_ali[$i]="-";
	#	    $template_ali[$i]="-";
	#	}
	#	elsif($first_aligned==0)
	#	{
	#	    $first_aligned=$i+1;
	#	}
	#	else
	#	{
	#	    $last_aligned=$i+1;
	#	}
	#    }
	#    my $target_str=join('',@target_ali[0..$#target_ali]);
	#    my $template_str=join('',@template_ali[0..$#template_ali]); 
	#    
	#    my $target_seq = new Bio::LocatableSeq(-moltype => "protein",
	#					   -seq => $target_str,
	#					   -id => "Target: ");
	#    my $template_seq = new Bio::LocatableSeq(-moltype => "protein",
	#					     -seq=> $template_str,
	#					     -id=> "$header");
	#    $aln->add_seq($target_seq);
	#    $aln->add_seq($template_seq); 
	#    return $aln;
	#}
    }
    #return $aln;
}
	

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in casp format
 Returns : 1 for success and 0 for error
 Args    : Bio::SimpleAlign object


=cut

sub write_aln {
    my ($self,@aln) = @_;
    my @seqs=();
    my ($rseq,$name,$count,$length,$seqsub);
   
  foreach my $aln (@aln) 
  {
      foreach $rseq ( $aln->each_seq() ) 
      {
	  $name = $rseq->id();  
	  push(@seqs,$rseq->seq());
      }
      #my @list=split(/\s+/,$name);
      #my $pdb_id=$list[1];
      $self->_print("PARENT $name\n");
      my @target_seq=split(//,$seqs[0]);
      my @template_seq=split(//,$seqs[1]);
      for(my $i=0;$i<scalar @target_seq;$i++)
      {
	  if($target_seq[$i] ne '-')
	  {
	      my $number=$i+1;
	      my $line="$target_seq[$i]\t$number\t$template_seq[$i]\t???\n";
	      $self->_print($line);
	  }
      }
      @seqs=();
      $self->_print("TER\n");
  }
    return 1;
}

1;
