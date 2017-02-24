# $Id: LocatableSeq.pm,v 1.10 2001/06/12 12:08:56 heikki Exp $
#
# BioPerl module for Bio::LocatableSeq
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::LocatableSeq - A Sequence object with start/end points on it

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION


    # a normal sequence object
    $locseq->seq();
    $locseq->id();

    # has start,end points
    $locseq->start();
    $locseq->end();

    # inheriets off RangeI, so range operations possible

=head1 FEEDBACK


=head2 Mailing Lists


User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists


The locatable sequence object was developed mainly because the 
SimpleAlign object requires this functionality, and in the rewrite
of the Sequence object we had to decide what to do with this.

It is, to be honest, not well integrated with the rest of bioperl, for
example, the trunc() function does not return a LocatableSeq object,
as some might have thought. There are all sorts of nasty gotcha's about
interactions between coordinate systems when these sort of objects are
used. 


=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...

package Bio::LocatableSeq;
use vars qw(@ISA);
use strict;

use Bio::PrimarySeq;
use Bio::RangeI;

@ISA = qw(Bio::PrimarySeq Bio::RangeI);

# new() is inherited from Bio::Root::RootI

# _initialize is where the heavy stuff will happen when new is called

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($start,$end,$strand) = $self->_rearrange( [qw(START END STRAND)],@args);
    
    defined $start && $self->start($start);
    defined $end   && $self->end($end);
    defined $strand && $self->strand($strand);
    
    return $self; # success - we hope!
}

=head2 start

 Title   : start
 Usage   : $obj->start($newval)
 Function: 
 Returns : value of start
 Args    : newvalue (optional)

=cut

sub start{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'start'} = $value;
    }
    return $self->{'start'};

}

=head2 end

 Title   : end
 Usage   : $obj->end($newval)
 Function: 
 Returns : value of end
 Args    : newvalue (optional)

=cut

sub end {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      my $string = $self->seq;
      if ($string and $self->start) {
	  my $s2 = $string;
	  $string =~ s/[.-]+//g;
	  my $len = CORE::length $string;
	  my $new_end = $self->start + $len - 1 ;
	  my $id = $self->id; 
	  $self->warn("In sequence $id residue count gives value $len.
Overriding value [$value] with value $new_end for Bio::LocatableSeq::end().")
	      and $value = $new_end if $new_end != $value;
      }
      $self->{'end'} = $value;
    }
    return $self->{'end'};

}

=head2 strand

 Title   : strand
 Usage   : $obj->strand($newval)
 Function: 
 Returns : value of strand
 Args    : newvalue (optional)

=cut

sub strand{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'strand'} = $value;
    }
    return $self->{'strand'};

}

=head2 get_nse

 Title   : get_nse
 Usage   :
 Function: read-only name of form id/start-end 
 Example :
 Returns : 
 Args    :

=cut

sub get_nse{
   my ($self,$char1,$char2) = @_;
  
   $char1 ||= "/";
   $char2 ||= "-";

   $self->throw("Attribute id not set") unless $self->id();
   $self->throw("Attribute start not set") unless $self->start();
   $self->throw("Attribute end not set") unless $self->end();

   return $self->id() . $char1 . $self->start . $char2 . $self->end ;

}


=head2 no_gaps

 Title   : no_gaps
 Usage   :$self->no_gaps('.')
 Function: 

           Gets number of gaps in the sequence. The count excludes
           leading or trailing gap characters.

           Valid bioperl sequence characters are [A-Za-z\-\.\*]. Of
           these, '.' and '-' are counted as gap characters unless an
           optional argument specifies one of them. 

 Returns : number of internal gaps in the sequnce. 
 Args    : a gap character (optional)

=cut

sub no_gaps {
    my ($self,$char) = @_;
    my ($seq, $count) = (undef, 0);

    # default gap characters
    $char ||= '-.';

    $self->warn("I hope you know what you are doing setting gap to [$char]")
	unless $char =~ /[-.]/;

    $seq = $self->seq;
    return 0 unless $seq; # empty sequence does not have gaps

    $seq =~ s/^([$char]+)//;
    $seq =~ s/([$char]+)$//;
    $count++ while $seq =~ /[$char]+/g;

    return $count;

}


=head2 column_from_residue_number

 Title   : column_from_residue_number
 Usage   : $col = $seq->column_from_residue_number($resnumber)
 Function:

           This function gives the position in the alignment
           (i.e. column number) of the given residue number in the
           sequence. For example, for the sequence

  	     Seq1/91-97 AC..DEF.GH

           column_from_residue_number(94) returns 5.

           An exception is thrown if the residue number would lie
           outside the length of the aligment
           (e.g. column_from_residue_number( "Seq2", 22 )

 Returns : A column number for the position of the
           given residue in the given sequence (1 = first column)
 Args    : A residue number in the whole sequence (not just that
           segment of it in the alignment)

=cut

sub column_from_residue_number {
    my ($self, $resnumber) = @_;

    $self->throw("Residue number has to be a positive integer, not [$resnumber]") 
	unless $resnumber =~ /^\d+$/ and $resnumber > 0;

    if ($resnumber >= $self->start() and $resnumber <= $self->end()) {
	my @residues = split //, $self->seq;
	my $count = $self->start();
	my $i;
	for ($i=0; $i < @residues; $i++) {
	    if ($residues[$i] ne '.' and $residues[$i] ne '-') {
		$count == $resnumber and last;
		$count++;
	    }		    
	}
	# $i now holds the index of the column. The actual column number is this index + 1
	
	return $i+1;
    }

    $self->throw("Could not find residue number $resnumber");

}


1;
