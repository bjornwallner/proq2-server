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

Bio::AlignIO::modhmm - modhmm sequence input/output stream

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


package Bio::AlignIO::modhmm;
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
    my $template_id;
    my ($logodds,$logpval);
    while(defined ($line = $self->_readline)) # This will read in the one line of file.
    {
	chomp($line);
	if($line=~/^Template\s+name/)
	{
	    my @temp=split(/\s+/,$line);
	    $parent=$temp[2];
	    $template_id=uc(substr($parent,3,6));
	}
	if($line=~/^Logodds\s+score/)
	{
	    my @temp=split(/\s+/,$line);
	    $logodds=$temp[3];
	}
	if($line=~/^log\(P-value\)/)
	{
	    my @temp=split(/\s+/,$line);
	    $logpval=$temp[2];
	}
	if($line=~/^Target\s+:/)
	{
	    my @temp=split(/\s+/,$line);
	    $target_ali=$temp[2];
	    $target_ali=~s/;//g;
	}
	if($line=~/^Template\:/)
	{
	    my @temp=split(/\s+/,$line);
	    $template_ali=$temp[1];
	    $template_ali=~s/;//g;
	}
	if($line=~/^END/)
	{
	    my $header="REMARK PARENT $template_id $parent¤";
	    $header.="REMARK SCORE $logpval $logodds¤";
	    $header.="REMARK METHOD modhmm¤";
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
