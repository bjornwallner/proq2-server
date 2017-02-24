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

Bio::AlignIO::casp_all - fasta sequence input/output stream

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
package Bio::AlignIO::casp_all2;

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
   # print "hello\n";
    my $self = shift;
    my $line;
    my ($parent);
    my @target_ali=();
    my @template_ali=();
    my $ali_length=0;
    my $aln =  Bio::SimpleAlign->new();
    my $jump_to_next;
    my $first_aligned=0;
    my $last_aligned=0;
    my $header="";
    my $template_start="";
    
    my $pdbfile=0;
    my $template_end="";
    my $target_start="";
    my $target_end="";
    my $target_str="";
    my $template_str="";
    my $model_nr="";
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
    while(defined ($line = $self->_readline)) # This will read in the one line of file.
    {
	#print "CASP: $line";
	if($line=~/^PARENT\s+(\w+)/)
	{
	    #print "CASP: $line";
	    $jump_to_next=0;
	    $parent=uc($1);
	    chomp($line);
	    $header.=$line."¤";
	}
	elsif($line=~/^ATOM/)
	{
	    $pdbfile=1;
	    $target_str.=$line;

	}
	elsif($line=~/^REMARK MODEL/)
	{
	    
	    chomp($line);
	    $header.=$line."¤";
	    my @temp=split(/\s+/,$line);
	    $model_nr=$temp[2];
	    
	}
	elsif($line=~/^REMARK/)
	{
	    chomp($line);
	    $header.=$line."¤";

	}
	elsif($line=~/^FILE/)
	{
	   # print "CASP2: $line\n";
	    chomp($line);
	    $header.=$line."¤";
	}
	elsif($line=~/^MODEL/)
	{
	    
	    chomp($line);
	    #$header="" if(not($header=~/FILE/));
	    $header.=$line."¤";
	    my @temp=split(/\s+/,$line);
	    $model_nr=$temp[1];
	    @target_ali=();
	    @template_ali=();
	}
	elsif($line=~/^(\w)\s+(\d+)\s+(\w)\s+(\d+)/) # Does the line begin with one character?
	{
	    #print "CASP: $line";
	    my $index=$2-1;
	    $target_ali[$index]=$1;
	    $template_ali[$index]=$3;
	    $ali_length++;
	}
	elsif($line=~/^TER/ && scalar @target_ali>0 && not($pdbfile)) # && $model_nr ne "")
	{
	    if($ali_length<5)
	    {
		$header.="REMARK TOO SHORT ALIGNMENT IGNORED\n";
	    }
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


	    $target_str=join('',@target_ali[0..$#target_ali]);
	    $template_str=join('',@template_ali[0..$#template_ali]); 
	    my $target_id="Target: $parent";
	    my $template_id="Template: $parent";
	    

	    #$template_start=1;
	    #$template_end=length($template_str);
	    #print "$template_end $template_str\n";
	    #$target_start=1;
	    #$target_end=length($target_str);


	    my $target_seq = new Bio::LocatableSeq('-seq' => $target_str,
						   '-id' => "Target: $parent $first_aligned $last_aligned");
						   						   
	    my $template_seq = new Bio::LocatableSeq('-seq'=> $template_str,
						     '-id'=> "$header");
	    $aln->add_seq($target_seq);
	    $aln->add_seq($template_seq); 
	    return $aln;
	}
	elsif($line=~/^TER/ && scalar @target_ali==0 && not($pdbfile))
	{
	    $header="";
	}
	elsif($line=~/^TER/ && $pdbfile)   # Special case for shgu lgscore2 could be 
	                                   # remove since we only take the model directly anyway
	{
	    my $pdb_code=lc(substr($parent,0,4));
	    my $pdbchain=" ";
	    $pdbchain=substr($parent,5,1) if(length($parent)>4);
	    $pdbfile="pdb$pdb_code.ent";
	    my $target_file="target_shgu.$model_nr.pdb";
	    print $target_file."\n";
	    open(*TARGET,">$target_file") || die "Cannot open $target_file\n";
	    print TARGET "$target_str";
	    close(*TARGET);
	    my $target_seq = new Bio::LocatableSeq('-seq' => "AAAAACCCCCC",
						   '-id' => "Target: $parent $first_aligned $last_aligned");
	    
	    my $template_seq = new Bio::LocatableSeq('-seq'=> "AAAAACCCCCC",
						     '-id'=> "$header");
	    $aln->add_seq($target_seq);
	    $aln->add_seq($template_seq); 
	    return $aln;
	}
	elsif($line=~/^TER/ && scalar @target_ali==0)
	{
	    $header.=$line;
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

1;
