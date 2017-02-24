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

Bio::AlignIO::casp - fasta sequence input/output stream

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


package Bio::AlignIO::casp;
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
    my ($parent);
    my @target_ali=();
    my @template_ali=();
    my $aln =  Bio::SimpleAlign->new();
    my $jump_to_next;
    my $first_aligned=0;
    my $last_aligned=0;
    while(defined ($line = $self->_readline)) # This will read in the one line of file.
    {
	if($line=~/^PARENT (\w+)/)
	{
	    $jump_to_next=0;
	    $parent=$1;
	}
	elsif($line=~/^REMARK LATE ERROR/ || ($line=~/^TER/ && scalar(@target_ali) == 0))
	{
	    $jump_to_next=1;
	}
	elsif($line=~/^(\w)\s+(\d+)\s+(\w)\s+(\d+)/ && not $jump_to_next) # Does the line begin with one character?
	{
	    my $index=$2-1;
	    $target_ali[$index]=$1;
	    $template_ali[$index]=$3;
	    
	}
	elsif($line=~/^TER/ && not($jump_to_next))
	{
	    for(my $i=0;$i<scalar @target_ali;$i++)
	    {
		if(not(defined($target_ali[$i])))
		{
		    $target_ali[$i]="-";
		    $template_ali[$i]="-";
		}
		elsif($first_aligned==0)
		{
		    $first_aligned=$i+1;
		}
		else
		{
		    $last_aligned=$i+1;
		}
	    }
	    my $target_str=join('',@target_ali[0..$#target_ali]);
	    my $template_str=join('',@template_ali[0..$#template_ali]); 
	    my $target_id="Target: $parent";
	    my $template_id="Template: $parent";
	    my $target_seq = new Bio::LocatableSeq('-seq' => $target_str,
						   '-id' => "Target: $parent $first_aligned $last_aligned");
	    my $template_seq = new Bio::LocatableSeq('-seq'=> $template_str,
						     '-id'=> "Template: $parent");
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
