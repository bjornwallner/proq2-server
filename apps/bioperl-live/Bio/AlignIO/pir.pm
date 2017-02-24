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

Bio::AlignIO::pir - fasta sequence input/output stream

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

#require "/home/ae/bjorn/modules/perl/bjornlib.pl";
package Bio::AlignIO::pir;
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
    my $parent="";;
    my $aln =  Bio::SimpleAlign->new();
    my $target_str="";
    my $template_str="";
    my $get_template_ali=0;
    my $get_target_ali=1;

    my %aa321=('ALA', 'A',
	'ARG', 'R',
	'ASN', 'N',
	'ASP', 'D',
	'CYS', 'C',
	'GLN', 'Q',
	'GLU', 'E',
	'GLY', 'G',
	'HIS', 'H',
	'ILE', 'I',
	'LEU', 'L',
	'LYS', 'K',
	'MET', 'M',
	'PHE', 'F',
	'PRO', 'P',
	'SER', 'S',
	'THR', 'T',
	'TRP', 'W',
	'TYR', 'Y',
	'VAL', 'V',
        'ASX', 'B',
        'GLX', 'Z',
        'XXX', 'A',
        'MSE', 'M',
        'FME', 'M',
        'PCA', 'E',
        '5HP', 'E',
        'SAC', 'S',
        'CCS', 'C');
#>P1;1ZIN
#structure:1ZIN: .  : .  : . :  : : : :
#PRTVAQAEALETMLADIGRKLD*
#
#>P1;target_name
#sequence:target_name: . : . : . : : : : :
#PKLRAQAEALEVLLRDVSVRPD*
    

    while(defined ($line = $self->_readline)) # This will read in the one line of file.
    {
	$line=~s/\r\n/\n/g;
#	print $line;
	chomp($line);
	if($line=~/^>P1;\s{0,1}$/)
	{
	    print "Assume target....\n";
	    $line.="target";
	}
	
	if($line=~/^>P1;\s{0,1}([\w\_\.]+)/)
	{
	    if($get_target_ali) 
	    {
		$target_str="";
		$template_str="";
		$parent=$1;
#		my @temp=split(/>P1;/,$line);

		#$parent=uc($temp[1]);
		#$parent=~s/\r//g;
		#print "TEST1 $parent $1\n";
		#print "REAL $1\n";
		if(length($parent) != 7) # && not($parent=~/^d/))
		{
		    $parent=uc($parent);
#		    $parent=~s/\r//g;
		    my $chain="";
		    my $id=lc(substr($parent,0,4));
		    my $chain_temp=substr($parent,4,1);
		    if(length($parent)>4)
		    {
			my $chain_temp=substr($parent,4,1);
			if(not($chain_temp=~/[A-Z]/))
			{
			    $chain_temp=substr($parent,5,1);
			}
			if($chain_temp=~/[A-Z]/)
			{
			    $chain=$chain_temp;
			    $parent=$id."_".$chain;
			}
			else
			{
			    $parent=$id;
			}
		    }
		}
		else #scop domain.
		{
		    #print "TEST $parent\n";
		}
		$get_template_ali=1;
		$get_target_ali=0;
	    }
	    else
	    {
		$get_target_ali=1;
		$get_template_ali=0;
	    }
	}
	elsif($get_template_ali)
	{
	    chomp($line);
	    $template_str.=uc($line);
	    $template_str=~s/\*//g;
	    #chop($template_str);
	    $template_str=~s/\r//g;
	    $template_str=~s/ //g;
	}
	elsif($get_target_ali)
	{
	    chomp($line);
	    $target_str.=uc($line);
	    #$target_str=~s/\*//g;
	    $target_str=~s/\r//g;
	    $target_str=~s/ //g;
	    if($target_str=~/\*/)
	    {
		$target_str=~s/\*//g;
		my $len_target=length($target_str);
		my $len_template=length($template_str);
		
		if($len_target>$len_template)
		{
		    $template_str.=_dashes($len_target-$len_template);
		}
		elsif($len_target<$len_template)
		{
		    $target_str.=_dashes($len_template-$len_target);
		}
		my $target_seq = new Bio::LocatableSeq('-seq' => $target_str,
						       '-id' => "Target: $parent");
		
		my $template_seq = new Bio::LocatableSeq('-seq'=> $template_str,
							 '-id'=> "PARENT $parent¤");
	#	print "$target_str\n$template_str\n$parent\n";
		$aln->add_seq($target_seq);
		$aln->add_seq($template_seq); 
		return $aln;
		
	    }
		
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
sub _get_pdb
{
    my $pdb_code=shift;
    my $PDBURL="ftp://ftp.rcsb.org/pub/pdb/data/structures/all/pdb/";
    my $OBSOLETE_PDBURL="ftp://ftp.rcsb.org/pub/pdb/data/structures/obsolete/pdb/";
    my $MODEL_PDBURL="ftp://ftp.rcsb.org/pub/pdb/data/structures/models/current/pdb/";
    my $file="pdb$pdb_code.ent";

    print "Trying ordinary directory\n";
    my $count=1;
    while(not(-e "$file.Z") && $count<=3)
    {
	print "Try $count: ncftpget $PDBURL$file.Z\n";
	`ncftpget $PDBURL$file.Z`;
	$count++;
    }
    if(-e "$file.Z")
    {
	`gunzip $file.Z`;
    }
    if(!-e "$file")
    {
	print "Trying theoretical models directory\n";
	my $subdir=substr($pdb_code,1,2);
	$subdir.="/";
	my $pdbfile_compressed=$file.".Z";
	my $pdbfile_dir=$subdir.$pdbfile_compressed;
	$count=1;
	while(not(-e "$file.Z") && $count<=3)
	{
	    print "Try $count: $MODEL_PDBURL$pdbfile_dir\n";
	    `ncftpget $MODEL_PDBURL$pdbfile_dir`;
	    $count++;
	}
	if(-e "$file.Z")
	{
	    `gunzip $file.Z`;
	}
    }
    if(!-e "$file")
    {
	print "Trying obsolete directory\n";
	my $subdir=substr($pdb_code,1,2);
	$subdir.="/";
	my $pdbfile_compressed=$file.".Z";
	my $pdbfile_obsolete=$subdir.$pdbfile_compressed;
	$count=1;
	while(not(-e "$file.Z") && $count<=3)
	{
	    print "Try $count: $OBSOLETE_PDBURL$pdbfile_obsolete\n";
	    `ncftpget $OBSOLETE_PDBURL$pdbfile_obsolete`;
	    $count++;
	}
	if(-e "$file.Z")
	{
	    `gunzip $file.Z`;
	}
	
    }
    if(-e "$file")
    {
	return 1;
    }
    else
    {
	return 0;
    }
}
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

1;
