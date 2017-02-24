# $Id: RichSeqI.pm,v 1.5 2001/02/26 11:35:01 lapp Exp $
#
# BioPerl module for Bio::Seq::RichSeqI
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::RichSeqI - RichSeq interface, mainly for database orientated sequences

=head1 SYNOPSIS

    @secondary   = $richseq->get_secondary_accessions;
    $division    = $richseq->division;
    $mol         = $richseq->molecule;
    @dates       = $richseq->get_dates; 
    $seq_version = $richseq->seq_version;  
    $pid         = $richseq->pid;

=head1 DESCRIPTION

This interface extends the Bio::SeqI interface to give additional functionality
to sequences with richer data sources, in particular from database sequences 
(EMBL, GenBank and Swissprot).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::RichSeqI;
use vars qw(@ISA);
use strict;
use Bio::SeqI;

@ISA = ('Bio::SeqI');


=head2 get_secondary_accessions

 Title   : get_secondary_accessions
 Usage   : 
 Function: Get the secondary accessions for a sequence.
 Example :
 Returns : an array of strings
 Args    : none


=cut

sub get_secondary_accessions{
   my ($self,@args) = @_;

   $self->throw("hit get_secondary_accessions in interface definition - error");

}


=head2 division

 Title   : division
 Usage   :
 Function: Get (and set, depending on the implementation) the divison for
           a sequence.

           Examples from GenBank are PLN (plants), PRI (primates), etc.
 Example :
 Returns : a string
 Args    :


=cut

sub division{
   my ($self,@args) = @_;

   $self->throw("hit division in interface definition - error");

}


=head2 molecule

 Title   : molecule
 Usage   :
 Function: Get (and set, depending on the implementation) the molecule
           type for the sequence.

           This is not necessarily the same as Bio::PrimarySeqI::moltype(),
           because it is databank-specific.
 Example :
 Returns : a string
 Args    :


=cut

sub molecule{
   my ($self,@args) = @_;

   $self->throw("hit molecule in interface definition - error");
}

=head2 pid

 Title   : pid
 Usage   :
 Function: Get (and set, depending on the implementation) the PID property
           for the sequence.
 Example :
 Returns : a string
 Args    :


=cut

sub pid {
   my ($self,@args) = @_;

   $self->throw("hit pid in interface definition - error");
}

=head2 get_dates

 Title   : get_dates
 Usage   :
 Function: Get (and set, depending on the implementation) the dates the
           databank entry specified for the sequence
 Example :
 Returns : an array of strings
 Args    :


=cut

sub get_dates{
   my ($self,@args) = @_;

   $self->throw("hit get_dates in interface definition - error");

}


=head2 seq_version

 Title   : seq_version
 Usage   :
 Function: Get (and set, depending on the implementation) the version string
           of the sequence.
 Example :
 Returns : a string
 Args    :


=cut

sub seq_version{
   my ($self,@args) = @_;

   $self->throw("hit seq_version in interface definition - error");

}

1;