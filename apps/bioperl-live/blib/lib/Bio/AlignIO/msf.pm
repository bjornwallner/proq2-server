# $Id: msf.pm,v 1.8 2001/06/14 12:32:53 heikki Exp $
#
# BioPerl module for Bio::AlignIO::msf

#	based on the Bio::SeqIO::msf module
#       by Ewan Birney <birney@sanger.ac.uk>
#       and Lincoln Stein  <lstein@cshl.org>
#
#       and the SimpleAlign.pm module of Ewan Birney
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself
# _history
# September 5, 2000
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::msf - msf sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::AlignIO class.

=head1 DESCRIPTION

This object can transform Bio::SimpleAlign objects to and from msf flat
file databases.

=head1 FEEDBACK

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHORS - Peter Schattner

Email: schattner@alum.mit.edu


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::msf;
use vars qw(@ISA %valid_type);
use strict;

use Bio::AlignIO;
use Bio::SeqIO::gcg; # for GCG_checksum()

@ISA = qw(Bio::AlignIO);

BEGIN {
    %valid_type = qw( dna N rna N protein P );
}

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream. Tries to read *all* MSF
          It reads all non whitespace characters in the alignment
          area. For MSFs with weird gaps (eg ~~~) map them by using
          $al->map_chars('~','-')
 Returns : SimpleAlign object
 Args    : NONE

=cut

sub next_aln {
    my $self = shift;
    my $entry;
    my (%hash,$name,$str,@names,$seqname,$start,$end,$count,$seq);

    my $aln =  Bio::SimpleAlign->new();


    while( $entry = $self->_readline) {
        $entry =~ /\/\// && last; # move to alignment section
        $entry =~ /Name:\s+(\S+)/ && do { $name = $1;
			       		$hash{$name} = ""; # blank line
			       		push(@names,$name); # we need it ordered!
			   		};
       # otherwise - skip
    }

   # alignment section

   while( $entry = $self->_readline) {
       $entry =~ /^\s*(\S+)\s+(.*)$/ && do {
	   $name = $1;
	   $str = $2;
	   if( ! exists $hash{$name} ) {
	       $self->throw("$name exists as an alignment line but not in the header. Not confident of what is going on!");
	   }
	   $str =~ s/\s//g;
	   $hash{$name} .= $str;
       };
   }

   return 0 if scalar @names < 1;

   # now got this as a name - sequence hash. Lets make some sequences!

   foreach $name ( @names ) {
       if( $name =~ /(\S+)\/(\d+)-(\d+)/ ) {
	   $seqname = $1;
	   $start = $2;
	   $end = $3;
       } else {
	   $seqname=$name;
	   $start = 1;
	   $str = $hash{$name};
	   $str =~ s/[^A-Za-z]//g;
	   $end = length($str);
       }

       $seq = new Bio::LocatableSeq('-seq'=>$hash{$name},
			   '-id'=>$seqname,
			   '-start'=>$start,
			   '-end'=>$end,
			   );

       $aln->add_seq($seq);


#  If $end <= 0, we have either reached the end of
#  file in <> or we have encountered some other error
#
#   if ($end <= 0) { undef $aln;}


   }

   return $aln;
}




=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in MSF format
           Sequence type of the alignment is determined by the first sequence.
 Returns : 1 for success and 0 for error
 Args    : Bio::SimpleAlign object


=cut

sub write_aln {
    my ($self,@aln) = @_;
    my $msftag;
    my $type;
    my $count = 0;
    my $maxname;
    my ($length,$date,$name,$seq,$miss,$pad,%hash,@arr,$tempcount,$index);
  foreach my $aln (@aln) {

    $date = localtime(time);
    $msftag = "MSF";
    $type = $valid_type{$aln->get_seq_by_pos(1)->moltype};
    $maxname = $aln->maxdisplayname_length();
    $length  = $aln->length();
    $name = $aln->id();
    if( !defined $name ) {
	$name = "Align";
    }


   $self->_print (sprintf("\n%s   MSF: %d  Type: %s  %s  Check: 00 ..\n\n", 
			  $name,  $aln->no_sequences, $type, $date));


      foreach $seq ( $aln->each_seq() ) {


	$name = $aln->displayname($seq->get_nse());
	$miss = $maxname - length ($name);
	$miss += 2;
	$pad  = " " x $miss;

	$self->_print (sprintf(" Name: %s%sLen:    %d  Check:  %d  Weight:  1.00\n",$name,$pad,length $seq->seq(), Bio::SeqIO::gcg->GCG_checksum($seq)));
	
	$hash{$name} = $seq->seq();
	push(@arr,$name);
    }
    	# ok - heavy handed, but there you go.
    	#
    	$self->_print ("\n//\n\n\n");

    	while( $count < $length ) {	
		# there is another block to go!
	   foreach $name  ( @arr ) {
	    	$self->_print (sprintf("%-20s  ",$name));
	
	    	$tempcount = $count;
	    	$index = 0;
	    	while( ($tempcount + 10 < $length) && ($index < 5)  ) {
		
			$self->_print (sprintf("%s ",substr($hash{$name},$tempcount,10)));
				
			$tempcount += 10;
			$index++;
	    	}	    	#
	    	# ok, could be the very last guy ;)
	    	#
	    	if( $index < 5) {
			# space to print!
			#
			$self->_print (sprintf("%s ",substr($hash{$name},$tempcount)));
			$tempcount += 10;
	    	}
	    	$self->_print ("\n");
  	   }
	    	$self->_print ("\n\n");
		$count = $tempcount;
    	}     			
      }
   return 1;
}

1;
