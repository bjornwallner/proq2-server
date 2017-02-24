#
# $Id: UTR3.pm,v 1.3 2001/06/19 23:02:16 mwilkinson Exp $
#

=head1 NAME

Bio::SeqFeature::UTR3

=head1 SYNOPSIS

  use Bio::SeqFeature::UTR3;
  my $Exon = Bio::SeqFeature::UTR3->new( # takes same args as SeqFeature::Generic
      -start => 10,
      -end => 100,
      -frame => '.',
      -primary => 'untranslated',
      -strand => '1',
      -source => 'cDNA_alignment',
      -score => '100',
      type => 'unknown');   # Descr. of the UTR3 type; *not* Cont.Vocab.


=head1 DESCRIPTION and ACKNOWLEDGEMENTS

Creates UTR3 type sequence features.  These are essentially SeaFeature::Generic
features, but report themselves as "Bio::SeqFeature::UTR3" when you query them
with a $Feature-E<gt>isa.

=head1 AUTHORS

Mark Wilkinson (mwilkinson@gene.pbi.nrc.ca)
Plant Biotechnology Institute, National Research Council of Canada.
Copyright (c) National Research Council of Canada, April, 2001.

=head1 DISCLAIMER

Anyone who intends to use and uses this software and code acknowledges and
agrees to the following: The National Research Council of Canada (herein "NRC")
disclaims any warranties, expressed, implied, or statutory, of any kind or
nature with respect to the software, including without limitation any warranty
or merchantability or fitness for a particular purpose.  NRC shall not be liable
in any event for any damages, whether direct or indirect,
consequential or incidental, arising from the use of the software.

=head1 CONTACT

Mark Wilkinson (mwilkinson@gene.pbi.nrc.ca)

=head1 METHODS

identical to SeqFeature::Generic except for:

=head2 type

  Usage:  $UTR3->type($type);
  Args:   optional string indicating new type
  Returns: current or newly set type


=cut


package Bio::SeqFeature::UTR3;

use strict;
use Carp;
use vars qw(@ISA $AUTOLOAD);
use Bio::SeqFeature::Generic;
@ISA = qw(Bio::SeqFeature::Generic);

{
	#class data	
	#___________________________________________________________
	#ATTRIBUTES (default)
    my %_attr_data = #     				DEFAULT    	ACCESSIBILITY
                  (
					type  		=> 	[undef,    'read/write'],  # should be something like 'init', 'intr', 'term', 'sngl' (GenScan exon types)					                  	
                    );
   #_____________________________________________________________
   #Class attribs and methods

    # Is a specified object attribute accessible in a given mode
    sub _accessible  {
    my ($self, $attr, $mode) = @_;
    $_attr_data{$attr}[1] =~ /$mode/
    }

    # Classwide default value for a specified object attribute
    sub _default_for {
    my ($self, $attr) = @_;
    $_attr_data{$attr}[0];
    }

    # List of names of all specified object attributes
    sub _standard_keys {
    keys %_attr_data;
	}
	

}

sub new {
	my ($caller, %args) = @_;
	
	#create a generic feature based on %args
	# is this a call to duplicate an object?
	my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;

    my $self = $caller->SUPER::new(%args);

    foreach my $attrname ( $self->_standard_keys ) {
    	if (exists $args{$attrname}) {
		$self->{$attrname} = $args{$attrname} }
    elsif ($caller_is_obj) {
		$self->{$attrname} = $caller->{$attrname} }
    else {
		$self->{$attrname} = $self->_default_for($attrname) }
    }   	

    return $self;

}


sub AUTOLOAD {
    no strict "refs";
    my ($self, $newval) = @_;

    $AUTOLOAD =~ /.*::(\w+)/;

    my $attr=$1;
    if ($self->_accessible($attr,'write')) {

	*{$AUTOLOAD} = sub {
	    if (defined $_[1]) { $_[0]->{$attr} = $_[1] }
	    return $_[0]->{$attr};
	};    ### end of created subroutine

###  this is called first time only
	if (defined $newval) {
	    $self->{$attr} = $newval
	}
	return $self->{$attr};

    } elsif ($self->_accessible($attr,'read')) {

	*{$AUTOLOAD} = sub {
	    return $_[0]->{$attr} }; ### end of created subroutine
	return $self->{$attr}  }


    # Must have been a mistake then...
    croak "No such method: $AUTOLOAD";
}

1;


