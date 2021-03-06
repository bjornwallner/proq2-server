# $Id: SeqAnalysisParserI.pm,v 1.6 2001/05/16 14:57:35 heikki Exp $
#
# BioPerl module for Bio::SeqAnalysisParserI
#
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>,
# and Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Jason Stajich, Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqAnalysisParserI - Sequence analysis output parser interface  

=head1 SYNOPSIS

    # get a SeqAnalysisParserI somehow, e.g. by
    my $parser = Bio::Factory::SeqAnalysisParserFactory->get_parser(
                            '-input' => 'inputfile', '-method' => 'genscan');
    while( my $feature = $parser->next_feature() ) {
	print "Feature from ", $feature->start, " to ", $feature->end, "\n";
    }

=head1 DESCRIPTION

SeqAnalysisParserI is a generic interface for describing sequence analysis
result parsers. Sequence analysis in this sense is a search for similarities
or the identification of features on the sequence, like a databank search or a
a gene prediction result.

The concept behind this interface is to have a generic interface in sequence
annotation pipelines (as used e.g. in high-throughput automated
sequence annotation). This interface enables plug-and-play for new analysis
methods and their corresponding parsers without the necessity for modifying
the core of the annotation pipeline. In this concept the annotation pipeline
has to rely on only a list of methods for which to process the results, and a
factory from which it can obtain the corresponding parser implementing this
interface.

See Bio::Factory::SeqAnalysisParserFactoryI and
Bio::Factory::SeqAnalysisParserFactory for interface and an implementation
of the corresponding factory.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                - General discussion
  http://bio.perl.org/MailList.html    - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Hilmar Lapp, Jason Stajich

Email Hilmar Lapp E<lt>hlapp@gmx.netE<gt>, Jason Stajich E<lt>jason@chg.mc.duke.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::SeqAnalysisParserI;
use strict;
use vars qw(@ISA);
use Bio::Root::RootI;
use Carp;
@ISA = qw(Bio::Root::RootI);

=head2 new

 Title   : new
 Remark  : Classes implementing this interface are expected to implement
           recognition of at least the following named parameters:
           -file     input file (alternative to -fh)
           -fh       input stream (alternative to -file)

=cut

sub new {
    my ($class,@args) = @_;
    if( $class eq 'Bio::SeqAnalysisParserI') {
	confess('Do not use this object directly it is an interface');
    }
    return $class->SUPER::new(@args);
}

=head2 next_feature

 Title   : next_feature
 Usage   : $seqfeature = $obj->next_feature();
 Function: Returns the next feature available in the analysis result, or
           undef if there are no more features.
 Example :
 Returns : A Bio::SeqFeatureI implementing object, or undef if there are no
           more features.
 Args    : none    

=cut

sub next_feature {
    my ($self);
    $self->_abstractDeath();
}

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  
  confess "Abstract method '$caller' defined in interface Bio::SeqAnalsysiParserI not implemented by pacakge $package. Not your fault - author of $package should be blamed!";
}

1;
