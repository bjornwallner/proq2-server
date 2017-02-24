# $Id: SimpleAlign.pm,v 1.28 2001/06/14 12:38:53 heikki Exp $
# BioPerl module for SimpleAlign
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code
#
#  History:
#	11/3/00 Added threshold feature to consensus and consensus_aa  - PS
#	May 2001 major rewrite - Heikki Lehvaslaiho

=head1 NAME

SimpleAlign - Multiple alignments held as a set of sequences

=head1 SYNOPSIS

  # use Bio::AlignIO to read in the alignment
  $str = Bio::AlignIO->new('-file' => 't/data/testaln.pfam');
  $aln = $str->next_aln();

  # some descriptors
  print $aln->length, "\n";
  print $aln->no_residues, "\n";
  print $aln->is_flush, "\n";
  print $aln->no_sequences, "\n";
  print $aln->percentage_identity, "\n";
  print $aln->consensus_string(50), "\n";

  # find the position in the alignment for a sequence location
  $pos = $aln->column_from_residue_number('1433_LYCES', 14); # = 6;

  # extract sequences and check values for the alignment column $pos
  foreach $seq ($aln->each_seq) {
      $res = $seq->subseq($pos, $pos);
      $count{$res}++;
  }
  foreach $res (keys %count) {
      printf "Res: %s  Count: %2d\n", $res, $count{$res};
  }


=head1 DESCRIPTION

SimpleAlign handles multiple alignments of sequences. It is very
permissive of types (it won\'t insist on things being all same length
etc): really it is a SequenceSet explicitly held in memory with a
whole series of built in manipulations and especially file format
systems for read/writing alignments.

SimpleAlign basically views an alignment as an immutable block of
text.  SimpleAlign *is not* the object to be using if you want to
perform complex alignment manipulations.  These functions are much
better done by UnivAln by Georg Fuellen.

However for lightweight display/formatting and minimal manipulation
(e.g. removing all-gaps columns) - this is the one to use.

SimpleAlign uses a subclass of L<Bio::PrimarySeq> class
L<Bio::LocatableSeq> to store its sequences. These are subsequences
with a start and end positions in the parent reference sequence.

Tricky concepts. SimpleAlign expects name,start,end to be 'unique' in
the alignment, and this is the key for the internal hashes.
(name,start,end is abbreviated nse in the code). However, in many
cases people don\'t want the name/start-end to be displayed: either
multiple names in an alignment or names specific to the alignment
(ROA1_HUMAN_1, ROA1_HUMAN_2 etc). These names are called
'displayname', and generally is what is used to print out the
alignment. They default to name/start-end.

The SimpleAlign Module came from Ewan Birney\'s Align module.

=head1 PROGRESS

SimpleAlign is being slowly converted to bioperl coding standards,
mainly by Ewan.

=over 3

=item Use Bio::Root::Object - done

=item Use proper exceptions - done

=item Use hashed constructor - not done!

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

   bioperl-l@bioperl.org             - General discussion
   http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via email
or the web:

    bioperl-bugs@bio.perl.org
    http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR

Ewan Birney, birney@sanger.ac.uk

=head1 SEE ALSO

L<Bio::LocatableSeq.pm>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SimpleAlign;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::LocatableSeq;         # uses Seq's as list

@ISA = qw(Bio::Root::RootI);

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  # we need to set up internal hashs first!

  $self->{'_seq'} = {};
  $self->{'_order'} = {};
  $self->{'start_end_lists'} = {};
  $self->{'_dis_name'} = {};
  $self->{'_id'} = 'NoName';

  # maybe we should automatically read in from args. Hmmm...

  return $self; # success - we hope!
}

=head1 Modifier methods

These methods modify the MSE by adding, removing or shuffling complete
sequences.

=head2 add_seq

 Title     : add_seq
 Usage     : $myalign->add_seq($newseq);
 Function  : Adds another sequence to the alignment. *Does not* align
             it - just adds it to the hashes.
 Returns   : nothing
 Args      : a L<Bio::LocatableSeq> object
             order (optional)

=cut

sub addSeq {
    my $self = shift;
    $self->warn(ref($self). "::addSeq - deprecated method. Use add_seq() instead.");
    $self->add_seq(@_);
}

sub add_seq {
    my $self = shift;
    my $seq  = shift;
    my $order = shift;
    my ($name,$id,$start,$end);

    if( !ref $seq || ! $seq->isa('Bio::LocatableSeq') ) {
	$self->throw("Unable to process non locatable sequences [", ref($seq), "]");
    }

    $id = $seq->id();
    $start = $seq->start();
    $end  = $seq->end();

    if( !defined $order ) {
	$order = keys %{$self->{'_seq'}};
    }

    $name = sprintf("%s/%d-%d",$id,$start,$end);

    if( $self->{'_seq'}->{$name} ) {
	$self->warn("Replacing one sequence [$name]\n");

    }
    else {
	#print STDERR "Assigning $name to $order\n";

	$self->{'_order'}->{$order} = $name;

	if (not exists( $self->{'_start_end_lists'}->{$id})) {
	    $self->{'_start_end_lists'}->{$id} = [];
	}
	push @{$self->{'_start_end_lists'}->{$id}}, $seq;
    }

    $self->{'_seq'}->{$name} = $seq;

}


=head2 remove_seq

 Title     : remove_seq
 Usage     : $aln->remove_seq($seq);
 Function  : Removes a single sequence from an alignment
 Returns   :
 Argument  : a Bio::LocatableSeq object

=cut

sub removeSeq {
    my $self = shift;
    $self->warn(ref($self). "::removeSeq - deprecated method. Use remove_seq() instead.");
    $self->remove_seq(@_);
}

sub remove_seq {
    my $self = shift;
    my $seq = shift;
    my ($name,$id,$start,$end);

    $self->throw("Need Bio::Locatable seq argument ")
	unless ref $seq eq 'Bio::LocatableSeq';

    $id = $seq->id();
    $start = $seq->start();
    $end  = $seq->end();
    $name = sprintf("%s/%d-%d",$id,$start,$end);

    if( !exists $self->{'_seq'}->{$name} ) {
	$self->throw("Sequence $name does not exist in the alignment to remove!");
    }

    delete $self->{'_seq'}->{$name};

    # we need to remove this seq from the start_end_lists hash

    if (exists $self->{'_start_end_lists'}->{$id}) {
	# we need to find the sequence in the array.
	
	my ($i, $found);;
	for ($i=0; $i < @{$self->{'_start_end_lists'}->{$id}}; $i++) {
	    if (${$self->{'_start_end_lists'}->{$id}}[$i] eq $seq) {
		$found = 1;
		last;
	    }
	}
	if ($found) {
	    splice @{$self->{'_start_end_lists'}->{$id}}, $i, 1;
	}
	else {
	    $self->throw("Could not find the sequence to remoce from the start-end list");
	}
    }
    else {
	$self->throw("There is no seq list for the name $id");
    }
    return 1;
    # we can't do anything about the order hash but that is ok
    # because each_seq will handle it
}


=head2 purge

 Title   : purge
 Usage   : $aln->purge(0.7);
 Function:

           Removes sequences above whatever %id.

           This function will grind on large alignments. Beware!
           (perhaps not ideally implemented)

 Example :
 Returns : An array of the removed sequences
 Args    :


=cut

sub purge{
  my ($self,$perc) = @_;
  my (@seqs,$seq,%removed,$i,$j,$count,@one,@two,$seq2,$k,$res,$ratio,@ret);

  @seqs = $self->each_seq();

#  foreach $seq ( @seqs ) {
#      printf("$seq %s %s\n",$seq->get_nse(),join(' ',$seq->dump()));
#  }

  for ($i=0;$i< @seqs;$i++ ) {
      $seq = $seqs[$i];

      # if it has already been removed, skip
      if( $removed{$seq->get_nse()} == 1 ) {
	  next;
      }

      @one = $seq->seq();
      for($j=$i+1;$j < @seqs;$j++) {
	  $seq2 = $seqs[$j];
	  if ( $removed{$seq2->get_nse()} == 1 ) {
	      next;
	  }
	  @two = $seq2->seq();
	  $count = 0;
	  $res = 0;
	  for($k=0;$k<@one;$k++) {
	      if( $one[$k] ne '.' && $one[$k] ne '-' && $one[$k] eq $two[$k]) {
		  $count++;
	      }
	      if( $one[$k] ne '.' && $one[$k] ne '-' && $two[$k] ne '.' && $two[$k] ne '-' ) {
		  $res++;
	      }
	  }
	  if( $res == 0 ) {
	      $ratio = 0;

	  } else {
	      $ratio = $count/$res;
	  }

	  if( $ratio > $perc ) {
	      $removed{$seq2->get_nse()} = 1;
	      $self->remove_seq($seq2);
	      push(@ret,$seq2);
	  } else {
	      # could put a comment here!
	  }
      }
  }
  
  return @ret;
}

=head2 sort_alphabetically

 Title     : sort_alphabetically
 Usage     : $ali->sort_alphabetically
 Function  : 

             Changes the order of the alignemnt to alphabetical on name 
             followed by numerical by number.

 Returns   : 
 Argument  : 

=cut

sub sort_alphabetically {
    my $self = shift;
    my ($seq,$nse,@arr,%hash,$count);

    foreach $seq ( $self->each_seq() ) {
	$nse = $seq->get_nse("-","-");
	$hash{$nse} = $seq;
    }

    $count = 0;

    %{$self->{'_order'}} = (); # reset the hash;

    foreach $nse ( sort _alpha_startend keys %hash) {
	$self->{'_order'}->{$count} = $nse;

	$count++;
    }

}

=head1 Sequence selection methods

Methods returning one or more sequences objects.

=head2 each_seq

 Title     : each_seq
 Usage     : foreach $seq ( $align->each_seq() ) 
 Function  : Gets an array of Seq objects from the alignment
 Returns   : an array
 Argument  : 

=cut

sub eachSeq {
    my $self = shift;
    $self->warn(ref($self). "::eachSeq - deprecated method. Use each_seq() instead.");
    $self->each_seq();
}

sub each_seq {
    my $self = shift;
    my (@arr,$order);

    foreach $order ( sort { $a <=> $b } keys %{$self->{'_order'}} ) {
	if( exists $self->{'_seq'}->{$self->{'_order'}->{$order}} ) {
	    push(@arr,$self->{'_seq'}->{$self->{'_order'}->{$order}});
	}
    }

    return @arr;
}


=head2 each_alphabetically

 Title     : each_alphabetically
 Usage     : foreach $seq ( $ali->each_alphabetically() )
 Function  :

             Returns an array of sequence object sorted alphabetically 
             by name and then by start point.
             Does not change the order of the alignment

 Returns   : 
 Argument  : 

=cut

sub each_alphabetically {
    my $self = shift;
    my ($seq,$nse,@arr,%hash,$count);

    foreach $seq ( $self->each_seq() ) {
	$nse = $seq->get_nse("-","-");
	$hash{$nse} = $seq;
    }

    foreach $nse ( sort _alpha_startend keys %hash) {
	push(@arr,$hash{$nse});
    }

    return @arr;

}

sub _alpha_startend {
    my ($aname,$astart,$bname,$bstart);
    ($aname,$astart) = split (/-/,$a);
    ($bname,$bstart) = split (/-/,$b);

    if( $aname eq $bname ) {
	return $astart <=> $bstart;
    }
    else { 
	return $aname cmp $bname;
    }

}

=head2 each_seq_with_id

 Title     : each_seq_with_id
 Usage     : foreach $seq ( $align->each_seq_with_id() ) 
 Function  : 

             Gets an array of Seq objects from the
             alignment, the contents being those sequences
             with the given name (there may be more than one)

 Returns   : an array
 Argument  : a seq name

=cut

sub eachSeqWithId {
    my $self = shift;
    $self->warn(ref($self). "::eachSeqWithId - deprecated method. Use each_seq_with_id() instead.");
    $self->each_seq_with_id(@_);
}

sub each_seq_with_id {
    my $self = shift;
    my $id = shift;

    $self->throw("Method each_seq_with_id needs a sequence name argument")
	unless defined $id;

    my (@arr, $seq);

    if (exists($self->{'_start_end_lists'}->{$id})) {
	@arr = @{$self->{'_start_end_lists'}->{$id}};
    }
    return @arr;
}

=head2 get_seq_by_pos

 Title     : get_seq_by_pos
 Usage     : $seq = $aln->get_seq_by_pos(3) # third sequence from the alignment
 Function  : 

             Gets a sequence based on its position in the alignment.
             Numbering starts from 1.  Sequence positions larger than
             no_sequences() will thow an error.

 Returns   : a Bio::LocatableSeq object
 Args      : positive integer for the sequence osition

=cut

sub get_seq_by_pos {
    my $self = shift;
    my ($pos) = @_;

    $self->throw("Sequence position has to be a positive integer, not [$pos]") 
	unless $pos =~ /^\d+$/ and $pos > 0;
    $self->throw("No sequence at position [$pos]") 
	unless $pos <= $self->no_sequences ;

    my $nse = $self->{'_order'}->{--$pos};
    return $self->{'_seq'}->{$nse};
}

=head1 Create new alignments

The result of these methods are horizontal or vertical subsets of the
current MSE.


=head2 select

 Title     : select
 Usage     : $aln2 = $aln->select(1, 3) # three first sequences
 Function  : 

             Creates a new alignment from a continuous subset of
             sequences.  Numbering starts from 1.  Sequence positions
             larger than no_sequences() will thow an error.

 Returns   : a Bio::SimpleAlign object
 Args      : positive integer for the first sequence
             positive integer for the last sequence to include (optional)

=cut

sub select {
    my $self = shift;
    my ($start, $end) = @_;

    $self->throw("Slice start has to be a positive integer, not [$start]") 
	unless $start =~ /^\d+$/ and $start > 0;
    $self->throw("Slice end has to be a positive integer, not [$end]") 
	unless $end  =~ /^\d+$/ and $end > 0;
    $self->throw("Slice $start [$start] has to be smaller than or equal to end [$end]") 
	unless $start <= $end;
    
    my $aln = new $self;
    foreach my $pos ($start .. $end) {
	$aln->add_seq($self->get_seq_by_pos($pos));
    }
    $aln->id($self->id);	
    return $aln;
}


=head2 slice

 Title     : slice
 Usage     : $aln2 = $aln->slice(20, 30)
 Function  : 

             Creates a slice from the alignment inclusive of start and
             end columns.  Sequences with no residues in the slice are
             excluded from the new alignment and a warning is printed.
             Slice beyond the length of the sequence does not do
             padding.

 Returns   : a Bio::SimpleAlign object
 Args      : positive integer for start column 
             positive integer for end column 

=cut

sub slice {
    my $self = shift;
    my ($start, $end) = @_;

    $self->throw("Slice start has to be a positive integer, not [$start]") 
	unless $start =~ /^\d+$/ and $start > 0;
    $self->throw("Slice end has to be a positive integer, not [$end]") 
	unless $end =~ /^\d+$/ and $end > 0;
    $self->throw("Slice $start [$start] has to be smaller than or equal to end [$end]") 
	unless $start <= $end;
    my $aln_length = $self->length;
    $self->throw("This alignment has only ". $self->length. 
		  " residues. Slice start [$start] is too bigger.") 
	 if $start > $self->length;

    my $aln = new $self;
    $aln->id($self->id);
    foreach my $seq ( $self->each_seq() ) {

	my $new_seq = new Bio::LocatableSeq (-id => $seq->id);

	# seq
	my $seq_end = $end;
	$seq_end = $seq->length if $end > $seq->length;
	my $slice_seq = $seq->subseq($start, $seq_end);
	$new_seq->seq( $slice_seq );
	
	# start
	if ($start > 1) {
	    my $pre_start_seq = $seq->subseq(1, $start - 1);
	    $pre_start_seq =~ s/\W//g; #print "$pre_start_seq\n";
	    $new_seq->start( $seq->start + CORE::length($pre_start_seq)  );
	} else {
	    $new_seq->start( $seq->start);
	}

	# end
	$slice_seq =~ s/\W//g;
	$new_seq->end( $new_seq->start + CORE::length($slice_seq) - 1 );

	if ($new_seq->start and $new_seq->end >= $new_seq->start) {
	    $aln->add_seq($new_seq);
	} else {
	    my $nse = $seq->get_nse();
	    $self->warn("Slice [$start-$end] of sequence [$nse] contains no residues. Sequence excluded from the new alignment.");
	}

    }

    return $aln;
}

=head1 Change sequences within the MSE

These methods affect characters in all sequences without changeing the
alignment.


=head2 map_chars

 Title     : map_chars
 Usage     : $ali->map_chars('\.','-')
 Function  : 

             Does a s/$arg1/$arg2/ on the sequences. Useful for gap
             characters

             Notice that the from (arg1) is interpretted as a regex,
             so be careful about quoting meta characters (eg
             $ali->map_chars('.','-') wont do what you want)

 Returns   : 
 Argument  : 'from' rexexp
             'to' string

=cut

sub map_chars {
    my $self = shift;
    my $from = shift;
    my $to   = shift;
    my ($seq,$temp);

    $self->throw("Need exactly two arguments") 
	unless defined $from and defined $to;

    foreach $seq ( $self->each_seq() ) {
	$temp = $seq->seq();
	$temp =~ s/$from/$to/g;
	$seq->seq($temp);
    }
    return 1;
}


=head2 uppercase

 Title     : uppercase()
 Usage     : $ali->uppercase()
 Function  : Sets all the sequences to uppercase
 Returns   : 
 Argument  : 

=cut

sub uppercase {
    my $self = shift;
    my $seq;
    my $temp;

    foreach $seq ( $self->each_seq() ) {
      $temp = $seq->seq();
      $temp =~ tr/[a-z]/[A-Z]/;

      $seq->seq($temp);
    }
    return 1;
}


=head2 match

 Title     : match()
 Usage     : $ali->match()
 Function  : 

             Goes through all columns and changes residues that are
             identical to residue in first sequence to match '.'
             character. Sets L<match_char>.

             USE WITH CARE: Most MSE formats do not support match
             characters in sequences, so this is mostly for output
             only. NEXUS format (L<Bio::AlignIO::nexus>) can handle
             it.

 Returns   : 1
 Argument  : a match character, optional, defaults to '.'

=cut

sub match {
    my ($self, $match) = @_;

    $match ||= '.';
    my ($matching_char) = $match;
    $matching_char = "\\$match" if $match =~ /[\^.$|()\[\]]/ ;  #'; 
    $self->map_chars($matching_char, '-');

    my @seqs = $self->each_seq();
    return 1 unless scalar @seqs > 1;

    my $refseq = shift @seqs ;
    my @refseq = split //, $refseq->seq;
    foreach my $seq ( @seqs ) {
	my @varseq = split //, $seq->seq();
	for ( my $i=0; $i < scalar @varseq; $i++) {
	    $varseq[$i] = $match if defined $refseq[$i] and 
		$refseq[$i] =~ /[A-Za-z\*]/ and $refseq[$i] eq $varseq[$i];
	}
	$seq->seq(join '', @varseq);
    }
    $self->match_char($match);
    return 1;
}


=head2 unmatch

 Title     : unmatch()
 Usage     : $ali->unmatch()
 Function  : 

             Undoes the effect of method L<match>. Unsets L<match_char>.

 Returns   : 1
 Argument  : a match character, optional, defaults to '.'

=cut

sub unmatch {
    my ($self, $match) = @_;

    $match ||= '.';

    my @seqs = $self->each_seq();
    return 1 unless scalar @seqs > 1;

    my $refseq = shift @seqs ;
    my @refseq = split //, $refseq->seq;
    foreach my $seq ( @seqs ) {
	my @varseq = split //, $seq->seq();
	for ( my $i=0; $i < scalar @varseq; $i++) {
	    $varseq[$i] = $refseq[$i] if defined $refseq[$i] and 
		$refseq[$i] =~ /[A-Za-z\*]/ and $varseq[$i] eq $match;
	}
	$seq->seq(join '', @varseq);
    }
    $self->match_char('');
    return 1;
}

=head1 MSE attibutes

Methods for setting and reading the MSE attributes. 

Note that the methods defining character semantics depend on the user
to set them sensibly.  They are needed only by certain input/output
methods. Unset them by setting to an empty string ('').

=head2 id

 Title     : id
 Usage     : $myalign->id("Ig")
 Function  : Gets/sets the id field of the alignment
 Returns   : An id string
 Argument  : An id string (optional)

=cut

sub id {
    my ($self, $name) = @_;

    if (defined( $name )) {
	$self->{'_id'} = $name;
    }
    
    return $self->{'_id'};
}

=head2 missing_char

 Title     : missing_char
 Usage     : $myalign->missing_char("?")
 Function  : Gets/sets the missing_char attribute of the alignment
             It is generally recommended to set it to 'n' or 'N' 
             for nucleotides and to 'X' for protein. 
 Returns   : An missing_char string,
 Argument  : An missing_char string (optional)

=cut

sub missing_char {
    my ($self, $char) = @_;

    if (defined $char ) {
	$self->throw("Single missing character, not [$char]!") if CORE::length($char) > 1;
	$self->{'_missing_char'} = $char;
    }
    
    return $self->{'_missing_char'};
}

=head2 match_char

 Title     : match_char
 Usage     : $myalign->match_char('.')
 Function  : Gets/sets the match_char attribute of the alignment
 Returns   : An match_char string,
 Argument  : An match_char string (optional)

=cut

sub match_char {
    my ($self, $char) = @_;

    if (defined $char ) {
	$self->throw("Single match character, not [$char]!") if CORE::length($char) > 1;
	$self->{'_match_char'} = $char;
    }

    return $self->{'_match_char'};
}

=head2 gap_char

 Title     : gap_char
 Usage     : $myalign->gap_char('-')
 Function  : Gets/sets the gap_char attribute of the alignment
 Returns   : An gap_char string, defaults to '-'
 Argument  : An gap_char string (optional)

=cut

sub gap_char {
    my ($self, $char) = @_;
    
    if (defined $char ) {
	$self->throw("Single gap character, not [$char]!") if CORE::length($char) > 1;
	$self->{'_gap_char'} = $char;
    } else {
	$self->{'_gap_char'} = '-';
    }
    
    return $self->{'_gap_char'};
}

=head1 Alignment descriptors

These read only methods describe the MSE in various ways.


=head2 consensus_string

 Title     : consensus_string
 Usage     : $str = $ali->consensus_string($threshold_percent)
 Function  : Makes a consensus
 Returns   : 
 Argument  : Optional treshold ranging from 0 to 100.  
                If consensus residue appears in fewer than threshold %
		of the sequences at a given location, consensus_string
		will return a "?" at that location rather than the
		consensus letter. (Default value = 0%)

=cut

sub consensus_string {
    my $self = shift;
    my $threshold = shift;
    my $len;
    my ($out,$count);

    $out = "";

    $len = $self->length;

    foreach $count ( 0 .. $len ) {
	$out .= $self->consensus_aa($count,$threshold);
    }
    return $out;
}


=head2 consensus_aa

 Title     : consensus_aa
 Usage     : $consensus_residue = $ali->consensus_aa($residue_number, $threshold_percent)
 Function  : Makes a consensus
 Returns   :
 Argument  : Optional treshold ranging from 0 to 100.  
                If consensus residue appears in fewer than threshold %
		of the sequences at a given location, consensus_string
		will return a "?" at that location rather than the
		consensus letter. (Default value = 0%)

=cut

sub consensus_aa {
    my $self = shift;
    my $point = shift;
    my $threshold_percent = shift || -1 ;
    my ($seq,%hash,$count,$letter,$key);

    foreach $seq ( $self->each_seq() ) {
	$letter = substr($seq->seq,$point,1);
	($letter =~ /\./) && next;
	# print "Looking at $letter\n";
	$hash{$letter}++;
    }
    my $number_of_sequences = $self->no_sequences();
    my $threshold = $number_of_sequences * $threshold_percent / 100. ;
    $count = -1;
    $letter = '?';

    foreach $key ( keys %hash ) {
	# print "Now at $key $hash{$key}\n";
	if( $hash{$key} > $count && $hash{$key} > $threshold) {
	    $letter = $key;
	    $count = $hash{$key};
	}
    }
    return $letter;
}


=head2 is_flush

 Title     : is_flush
 Usage     : if( $ali->is_flush() )
           : 
           :
 Function  : Tells you whether the alignment 
           : is flush, ie all of the same length
           : 
           :
 Returns   : 1 or 0
 Argument  : 

=cut

sub is_flush {
    my $self = shift;
    my $seq;
    my $length = (-1);
    my $temp;
    
    foreach $seq ( $self->each_seq() ) {
	if( $length == (-1) ) {
	    $length = CORE::length($seq->seq());
	    next;
	}

	$temp = CORE::length($seq->seq());

	if( $temp != $length ) {
	    return 0;
	}
    }

    return 1;
}


=head2 length

 Title     : length()
 Usage     : $len = $ali->length() 
 Function  : Returns the maximum length of the alignment.
             To be sure the alignment is a block, use is_flush
 Returns   : 
 Argument  : 

=cut

sub length_aln {
    my $self = shift;
    $self->warn(ref($self). "::length_aln - deprecated method. Use length() instead.");
    $self->length(@_);
}

sub length {
    my $self = shift;
    my $seq;
    my $length = (-1);
    my ($temp,$len);

    foreach $seq ( $self->each_seq() ) {
	$temp = CORE::length($seq->seq());
	if( $temp > $length ) {
	    $length = $temp;
	}
    }

    return $length;
}


=head2 maxdisplayname_length

 Title     : maxdisplayname_length
 Usage     : $ali->maxdisplayname_length()
 Function  : 

             Gets the maximum length of the displayname in the
             alignment. Used in writing out various MSE formats.

 Returns   : integer
 Argument  : 

=cut

sub maxname_length {
    my $self = shift;
    $self->warn(ref($self). "::maxname_length - deprecated method. Use maxname_length() instead.");
    $self->maxdisplayname_length();
}

sub maxnse_length {
    my $self = shift;
    $self->warn(ref($self). "::maxnse_length - deprecated method. Use maxnse_length() instead.");
    $self->maxdisplayname_length();
}

sub maxdisplayname_length {
    my $self = shift;
    my $maxname = (-1);
    my ($seq,$len);

    foreach $seq ( $self->each_seq() ) {
	$len = CORE::length $self->displayname($seq->get_nse());

	if( $len > $maxname ) {
	    $maxname = $len;
	}
    }

    return $maxname;
}

=head2 no_residues

 Title     : no_residues
 Usage     : $no = $ali->no_residues
 Function  : number of residues in total in the alignment
 Returns   : integer
 Argument  : 

=cut

sub no_residues {
    my $self = shift;
    my $count = 0;

    foreach my $seq ($self->each_seq) {
	my $str = $seq->seq();

	$count += ($str =~ s/[^A-Za-z]//g);
    }

    return $count;
}

=head2 no_sequences

 Title     : no_sequences
 Usage     : $depth = $ali->no_sequences
 Function  : number of sequence in the sequence alignment
 Returns   : integer
 Argument  : 

=cut

sub no_sequences {
    my $self = shift;

    return scalar($self->each_seq);
}

=head2 percentage_identity

 Title   : percentage_identity
 Usage   : $id = $align->percentage_identity
 Function: The function uses a fast method to calculate the average 
           percentage identity of the alignment
 Returns : The average percentage identity of the alignment
 Args    : None

=cut

sub percentage_identity{
   my ($self,@args) = @_;

   my @alphabet = ('A','B','C','D','E','F','G','H','I','J','K','L','M',
                   'N','O','P','Q','R','S','T','U','V','W','X','Y','Z');

   my ($len, $total, $subtotal, $divisor, $subdivisor, @seqs, @countHashes);

   if (! $self->is_flush()) {
       $self->throw("All sequences in the alignment must be the same length");
   }

   @seqs = $self->each_seq();
   $len = $self->length();

   # load the each hash with correct keys for existence checks
   for( my $index=0; $index < $len; $index++) {
       foreach my $letter (@alphabet) {
	   $countHashes[$index]->{$letter} = 0;
       }
   }
   foreach my $seq (@seqs)  {
       my @seqChars = split //, $seq->seq();
       for( my $column=0; $column < @seqChars; $column++ ) {
	   my $char = uc($seqChars[$column]);
	   if (exists $countHashes[$column]->{$char}) {
	       $countHashes[$column]->{$char}++;
	   }
       }
   }

   $total = 0;
   $divisor = 0;
   for(my $column =0; $column < $len; $column++) {
       my %hash = %{$countHashes[$column]};
       $subdivisor = 0;
       foreach my $res (keys %hash) {
	   $total += $hash{$res}*($hash{$res} - 1);
	   $subdivisor += $hash{$res};
       }
       $divisor += $subdivisor * ($subdivisor - 1);
   }
   return ($total / $divisor )*100.0;
}

=head1 Alignment positions

A method to map a sequence position into an alignment column.

=head2 column_from_residue_number

 Title   : column_from_residue_number
 Usage   : $col = $ali->column_from_residue_number( $seqname, $resnumber)
 Function:

           This function gives the position in the alignment
           (i.e. column number) of the given residue number in the
           sequence with the given name. For example, for the
           alignment

  	     Seq1/91-97 AC..DEF.GH
  	     Seq2/24-30 ACGG.RTY..
  	     Seq3/43-51 AC.DDEFGHI

           column_from_residue_number( "Seq1", 94 ) returns 5.
           column_from_residue_number( "Seq2", 25 ) returns 2.
           column_from_residue_number( "Seq3", 50 ) returns 9.

           An exception is thrown if the residue number would lie
           outside the length of the aligment
           (e.g. column_from_residue_number( "Seq2", 22 )

	  Note: If the the parent sequence is represented by more than
	  one alignment sequence and the residue number is present in
	  them, this method finds only the first one.

 Returns : A column number for the position in the alignment of the
           given residue in the given sequence (1 = first column)
 Args    : A sequence id/name (not a name/start-end)
           A residue number in the whole sequence (not just that
           segment of it in the alignment)

=cut

sub column_from_residue_number {
    my ($self, $name, $resnumber) = @_;

    $self->throw("No sequence with name [$name]") unless $self->{'_start_end_lists'}->{$name};
    $self->throw("Second argument residue number missing") unless $resnumber;

    foreach my $seq ($self->each_seq_with_id($name)) {
	my $col;
	eval {
	    $col = $seq->column_from_residue_number($resnumber);
	};
	next if $@;		
	return $col;
    }

    $self->throw("Could not find a sequence segment in $name containing residue number $resnumber");

}

=head1 Sequence names

Methods to manipulate the display name. The default name based on the
sequence id and subsequence positions can be overridden in various
ways.

=head2 displayname

 Title     : displayname
 Usage     : $myalign->displayname("Ig", "IgA")
 Function  : Gets/sets the display name of a sequence in the alignment
           :
 Returns   : A display name string
 Argument  : name of the sequence
             displayname of the sequence (optional)

=cut

sub get_displayname {
    my $self = shift;
    $self->warn(ref($self). "::get_displayname - deprecated method. Use displayname() instead.");
    $self->displayname(@_);
}

sub set_displayname {
    my $self = shift;
    $self->warn(ref($self). "::set_displayname - deprecated method. Use displayname() instead.");
    $self->displayname(@_);
}

sub displayname {
    my ($self, $name, $disname) = @_;

    $self->throw("No sequence with name [$name]") unless $self->{'_seq'}->{$name};

    if(  $disname and  $name) {
	$self->{'_dis_name'}->{$name} = $disname;
	return $disname;
    }
    elsif( defined $self->{'_dis_name'}->{$name} ) {
	return  $self->{'_dis_name'}->{$name};
    } else {
	return $name;
    }
}

=head2 set_displayname_count

 Title     : set_displayname_count
 Usage     : $ali->set_displayname_count
 Function  : 

             Sets the names to be name_# where # is the number of
             times this name has been used.

 Returns   : 
 Argument  : 

=cut

sub set_displayname_count {
    my $self= shift;
    my (@arr,$name,$seq,$count,$temp,$nse);

    foreach $seq ( $self->each_alphabetically() ) {
	$nse = $seq->get_nse();

	#name will be set when this is the second
	#time (or greater) is has been seen

	if( defined $name and $name eq ($seq->id()) ) {
	    $temp = sprintf("%s_%s",$name,$count);
	    $self->displayname($nse,$temp);
	    $count++;
	} else {
	    $count = 1;
	    $name = $seq->id();
	    $temp = sprintf("%s_%s",$name,$count);
	    $self->displayname($nse,$temp);
	    $count++;
	}
    }
    return 1;
}

=head2 set_displayname_flat

 Title     : set_displayname_flat
 Usage     : $ali->set_displayname_flat()
 Function  : Makes all the sequences be displayed as just their name,
             not name/start-end
 Returns   : 1
 Argument  : 

=cut

sub set_displayname_flat {
    my $self = shift;
    my ($nse,$seq);

    foreach $seq ( $self->each_seq() ) {
	$nse = $seq->get_nse();
	$self->displayname($nse,$seq->id());
    }
    return 1;
}

=head2 set_displayname_normal

 Title     : set_displayname_normal
 Usage     : $ali->set_displayname_normal() 
 Function  : Makes all the sequences be displayed as name/start-end
 Returns   : 
 Argument  : 

=cut

sub set_displayname_normal {
    my $self = shift;
    my ($nse,$seq);

    foreach $seq ( $self->each_seq() ) {
	$nse = $seq->get_nse();
	$self->displayname($nse,$nse);
    }
    return 1;
}

1;
