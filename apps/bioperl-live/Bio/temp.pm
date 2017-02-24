=Head NAME

Bio::Pdb - Object for handling pdb-files

=head1 VERSION

=head1 SYNOPSIS

=head1 DESCRIPTION

=head2 Overview

=head1 AUTHOR

Björn Wallner
bjorn@sbc.su.se

=head1 COPYRIGHT

Copyright (c)  2001, Björn Wallner. All Rights reserved.
This module is free software. It may be used, redistributed
and/or modified under the same terms as Perl itself.

=cut

package Bio::Pdb;
use lib '/afs/nada.kth.se/home/f96/f96-bwa/bjorn/bioperl/bioperl-live';
use strict;
use Bio::Seq;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::Tools::pSW;

#use lib '/afs/nada.kth.se/home/f96/f96-bwa/bjorn/bioperl/bioperl-live/';

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


my %aa123=('A','ALA',
	   'R','ARG', 
	   'N','ASN', 
	   'D','ASP',
	   'C','CYS',
	   'Q','GLN', 
	   'E','GLU', 
	   'G','GLY', 
	   'H','HIS',
	   'I','ILE',
	   'L','LEU', 
	   'K','LYS', 
	   'M','MET', 
	   'F','PHE',
	   'P','PRO',
	   'S','SER', 
	   'T','THR', 
	   'W','TRP',
	   'Y','TYR', 
	   'V','VAL',
           'B','ASX',
           'Z','GLX');
=head2 new

   Title    : new
   Usage    : $pdb    = Bio::Pdb->new( -file => 'pdbXXXX.ent',
                                       -chain => 'A',
                                       -pdbcode => 'XXXX');

   Function : Returns a new pdb object for
              the pdb file defined by the file string or from the pdbcode
   Returns  : a new Bio::Pdb object

=cut

sub new 
{
    my $FTPURL="ftp://ftp.ebi.ac.uk/pub/databases/pdb/all_entries/compressed_files/";
    my($caller,@args) = @_;
    my $class = ref($caller) || $caller;
    my $self = {};
    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param; # lowercase keys
 
    if(defined($param{-file}))
    {
	if(defined $param{-chain})
	{
	    $self=_read_pdb($param{-file},uc($param{-chain}));
	}
	else
	{
	    $self=_read_pdb($param{-file}," ");
	}
    }
    elsif(defined($param{-pdbcode}))
    {
	my $file="pdb$param{-pdbcode}.ent";
	`ncftpget $FTPURL$file`;
	if(defined $param{-chain})
	{
	    $self=_read_pdb($file,uc($param{-chain}));
	}
	else
	{
	    $self=_read_pdb($file," ");
	}
	`rm $file`;
	
    }
    
    if($self)
    {
	bless ($self,$class);
	return $self;
    }
    else
    {
	return 0;
    }
}
=head2 _read_pdb

  Title   : _read_pdb
  Usage   : Internal function used by the new function.
  Returns : A a pdb object to be blessed.

=cut

sub _read_pdb 
{
    my ($file,$chain)=@_;
    my $no_seqres=1;   # Deals with case when the pdb SEQRES is missing
    my @list=split(/\//,$file);
    my $name=$list[length(@list)-1];
    my @pdbfield=();
    my @pdbatomno=();
    my @pdbatomtype=();
    my @pdbatom_alt_loc=();
    my @pdbres=();
    my @pdbresno=();
    my @pdbinsertion_code=();
    my @x=();
    my @y=();
    my @z=();
    my @resbegin_temp=();
    my @resCA_temp=();
    my @resbegin=();
    my @resCA=();
    my $seqres="";
    my $seqca="";
    my $SeqRes;
    my $SeqCA;
    my $SeqRes_ali;
    my $SeqCA_ali;
    my $atomcount=0;
    my $resname="";
    my $old_resname="undef";
    my $last_added="";
    my @PDBHEADER=('^HEADER','^COMPND','^SOURCE','^AUTHOR','^TITLE',
		   ,'^REVDAT','REMARK', '^JRNL','^HET ', '^HETNAM','^MODRES',
		   ,'^HELIX','^SHEET','^LOOP','^TURN', 'MODEL');    # Specifies which parts of the pdb header to take.
    my @pdbhead=();
    my $seq_begin="";   # A sequence for the residue of first atom in a residue,
                        # for doing an alignment with CA_res to determine which residues to
                        # print out.
    #my $success=1;
    #open(FILE,"$file") or $success=0;
    #my $done=0;
    #if($success)
    #{
    my $chain_found=0;
    my $first_chain="undef";
    if($chain eq " ") # if the specified chain is " ", it checks that, that chain exists and if it doesn't it takes the first chain.
    {
	open (FILE,"$file")|| die "Cannot open $file.\n";
	while(<FILE>)     # if the specified chain is " " Checking whether the specified chain does exist if it doesn't take the first
	{
	    if(/^ATOM/)
	    {
		my $current_chain=substr($_,21,1);
		#print $current_chain."\n";
		$first_chain=$current_chain if($first_chain eq "undef");
		if($current_chain eq $chain)
		{
		    $chain_found=1;
		    last;
		}
	    }

	}
	if($first_chain eq "undef")
	{
	    warn "The chain \"$chain\" does not exist in the pdb file: $file\nNo object created...Bla..\n";
	    return 0;
	}
	elsif(not($chain_found))
	{
	    warn "Specified chain: \" \" not found taking the first chain $first_chain instead.\n";
	    $chain=$first_chain;
	}
	
    }
    open(FILE,"$file") || die "Cannot open $file.\n";
    my $reached_TER=0;
    while(<FILE>)
    {
	chomp;
	foreach my $key (@PDBHEADER)
	{
	    if(/$key/)
	    {
		push(@pdbhead,$_);
		last;
	    }
	}
	
	if(/^SEQRES/)
	{
	    my $current_chain=substr($_,11,1);
	    if($current_chain eq $chain)
	    {
		$no_seqres=0;
		my $atoms=substr($_, 19, 52);
		my @list=split(/\s+/,$atoms);
		for(my $i=0;$i<=$#list;$i++)
		{
		    $seqres.=$aa321{$list[$i]} if(defined($aa321{$list[$i]}));
		}	
		
	    }
	}
	if(not($reached_TER) && (/^ATOM/) || (/^HETATM/ && (/MSE/ || /PCA/ || /5HP/ || /SAC/ || /CCS/ || /FME/)))
	{
	    my $current_chain=substr($_,21,1);
	    if($current_chain eq $chain)
	    {
		my $field = substr($_,0,6);
		$field=~s/ //g;
		my $atomno=substr($_, 7, 4);
		my $atomtype=substr($_, 13, 3);
		my $alt_loc=substr($_,16,1);
		my $res=substr($_,17, 3);
		my $resno=substr($_, 22, 4);
		my $insertion_code=substr($_,26,1);
	      
		$atomno=~s/ //g;
		$res=~s/ //g;
		if(defined($aa321{$res}))
		{
		    $resno=~s/ //g;
		    $pdbfield[$atomcount]=$field;
		    $pdbatomno[$atomcount]=$atomno;
		    $pdbatomtype[$atomcount]=$atomtype;
		    $pdbatom_alt_loc[$atomcount]=$alt_loc;
		    $pdbres[$atomcount]=$res;
		    $pdbresno[$atomcount]=$resno;
		    $pdbinsertion_code[$atomcount]=$insertion_code;
		    $x[$atomcount]=substr($_,30,8);
		    $y[$atomcount]=substr($_,38,8);
		    #print "$res $chain $resno $y[$atomcount]\n";
		    $z[$atomcount]=substr($_,46,8);
		    $x[$atomcount]=~s/ //g;
		    $y[$atomcount]=~s/ //g;
		    $z[$atomcount]=~s/ //g;
		    #print $_."\n" if($atomtype eq "CA ");
		    $resname="$resno$insertion_code";

		    #THIS TRANSLATES ALL HETATM TO ATOM AND TO THE CORRESPONDING RESIDUE
		    if($field=~/HETATM/)
		    {
			$pdbres[$atomcount]=$aa123{$aa321{$res}};
			$pdbfield[$atomcount]="ATOM";
		    }
		    if($resname ne $old_resname || $atomcount == 0)
			#if($atomtype eq "N  " || $atomcount == 0)
		    {
			#print "$atomtype\n";
			push(@resbegin_temp,$atomcount);
			$seq_begin.=$aa321{$res}; #store the "FIRST ATOM IN RES" sequence for which coordinates exists.
			#print "$resno $atomcount $atomtype";
		    }
		    if($atomtype eq "CA " && $last_added ne $resname)
			#if($atomtype eq "CA ")
		    {
			if(scalar @resCA_temp != 0)
			{
			    #print $_."\n";
			    my $last_index=$#resCA_temp;
			    my $dist=distance($x[$atomcount],$y[$atomcount],$z[$atomcount],$x[$resCA_temp[$#resCA_temp]],$y[$resCA_temp[$#resCA_temp]],$z[$resCA_temp[$#resCA_temp]]);
			    if($dist>8 || $dist < 2)
			    {
				warn "Distance between CA $last_added $resname is $dist Å.\n";
				#print "($x[$atomcount],$y[$atomcount],$z[$atomcount],$x[$resCA_temp[$#resCA_temp]],$y[$resCA_temp[$#resCA_temp]],$z[$resCA_temp[$#resCA_temp]])\n";
				#print $_."\n";
			    }
			}
			#print "$x[$atomcount],$y[$atomcount],$z[$atomcount],$x[$resCA_temp[$#resCA]],$y[$resCA_temp[$#resCA]],$z[$resCA_temp[$#resCA]]"."\n";
			$seqca.=$aa321{$res};  #store the CA sequence for which coordinates exists.
			push(@resCA_temp,$atomcount);
			#print " $resno $atomcount $atomtype".scalar @resCA_temp." ".scalar @resbegin_temp."\n";
			$last_added=$resname;
		    }
		    $atomcount++;
		    $old_resname=$resname;
		}
	    }
	}
	$reached_TER=1 if((/^TER/ || /^ENDMDL/) && $atomcount>0);   # In order to avoid reading HETATM after termination of residue coordinates.
    }
    
    if($atomcount == 0)
    {
	warn "The chain \"$chain\" does not exist in the pdb file: $file\nNo object created.\n";
	return 0;
    }
    else
    {
   			
	# Aligment between CA res and "first atom res", to filter out the
	# cases when CA is missing, when CA is missing but some other atom is present that atom
	# is put into the CA position.
		
	my($Seq_begin_ali2,$SeqCA_ali2)=_align($seq_begin,$seqca);
	#print "$Seq_begin_ali2\n\n$SeqCA_ali2\n\n";
	my $len1 = scalar @resbegin_temp;
	my $len2 = scalar @resCA_temp;
	if($len1 != $len2)
	{
	    warn "CA is probably missing in some residue $file $len1 $len2";
	    #print "$SeqRes_ali\n\n$SeqCA_ali\n\n";
	}
	    
	my @ali=split(//,$SeqCA_ali2);  # Aligment between CA res and "first atom res", to filter out the
	                                # cases when CA is missing.
	my @ali_begin=split(//,$Seq_begin_ali2);
	my $CA_index=0;
	my $resbegin_index=0;
	for(my $i=0;$i<scalar @ali;$i++)
	{
	    if($ali[$i] eq "-" && $ali_begin[$i] ne "-")  #if CA not aligned with "first atom" in resdiue change CA pointer to $resbegin
	    {
		    #changed 2002-01-29 /BW
		$resbegin[$i]=$resbegin_temp[$resbegin_index];
		$resCA[$i]=$resbegin_temp[$resbegin_index];  #"-";
		$ali[$i]=$ali_begin[$i];  #OBSERVE this adds another residue to seqca which in fact is not CA but some other atom belonging to that res.
		$resbegin_index++;
		
		#$resbegin[$i]="-";
		#$resCA[$i]="-";
		#print "$resbegin[$i] $resCA[$i]\n";
	    }
	    else
	    {
		$resbegin[$i]=$resbegin_temp[$resbegin_index];
		$resCA[$i]=$resCA_temp[$CA_index];
		#print "$resbegin[$i] $resCA[$i]\n";
		#my $temp=$aa321{$pdbres[$resbegin_temp[$j]]};
		$CA_index++;
		$resbegin_index++;
	    }
	}

	$seqca=join('',@ali);
	#print $seqca."\n";
	$SeqCA = Bio::Seq->new( -seq => $seqca,
				-id  => 'SEQCA',);	
	
	if($seqres)   #Produce Alignment
	{
	    $SeqRes = Bio::Seq->new( -seq => $seqres,
				 -id  => 'SEQRES',);
	    ($SeqRes_ali,$SeqCA_ali)=_align($seqres,$seqca);
	}
	else
	{
	    $SeqRes=0;
	    $SeqRes_ali=0;
	    $SeqCA_ali=$seqca;
	}

	
    }
    my $self={name => $name,
	      header => [@pdbhead],
	      chain => $chain,
	      seqres => $SeqRes,
	      seqCA => $SeqCA,
	      seqres_ali => $SeqRes_ali,
	      seqCA_ali => $SeqCA_ali,
	      field => [@pdbfield],
	      atomno => [@pdbatomno],
	      atomtype => [@pdbatomtype],
	      atom_alt_loc=> [@pdbatom_alt_loc],
	      res => [@pdbres],
	      resno => [@pdbresno],
	      insertion_code => [@pdbinsertion_code],
	      x => [@x],
	      y => [@y],
	      z => [@z],
	      resbegin => [@resbegin],   #NO:pointer to begin of a residue if == "-" CA is missing for that residues
	      resCA => [@resCA]};        #NO:pointer to CA of a residue if == "-" CA is missing for that residues
	      #table321 => {%aa321}};

    return $self;
}

sub get_chain
{
    my $self=shift;
    return $self->{chain};

}

sub get_ali
{
    my $self=shift;
    my $coord;
    my @list=$self->{z};
}

sub number_of_res
{
    my $self=shift;
    return length($self->{seqCA}->seq());
}


sub CAres
{
    my $self=shift;
    return $self->{seqCA};
}
#sub CAresname
#{
#    my $self=shift;
 #   my @resno=();
 #   if(length($self->{seqCA_ali}) > 1)
  #  {
#	my $len=length($self->{seqCA_ali});
#	my @ali=split('',$self->{seqCA_ali});

# Function which returns a vector with all 


sub get_seqres
{
    my $self=shift;
    return $self->{seqres}->seq();
}

sub get_seqres_CA_ali
{
    my $self=shift;
    return ($self->{seqres_ali},$self->{seqCA_ali});
    
}



=head2 CAresname

    Title   : CAresname
    Usage   : Helper function for fix_numbering

    Function: 
    Example : @array=$pdb_obj->CAresname()
    Returns : An array with all
    Args    : None besides the object of course.

=cut

sub CAresname
{
    my $self=shift;
    my @resno=();
    #if(scalar @{$self->{resCA}} > 1)
    #{
	#my $len=length($self->{resCA});
    my @ali=@{$self->{resCA}};
    for(my $i=0;$i<=$#ali;$i++)
    {
	if($ali[$i] ne "-")
	{
	    push(@resno,$self->{resno}[$self->{resCA}[$i]].$self->{insertion_code}[$self->{resCA}[$i]]);
	    #print $self->{res}[$self->{resCA}[$i]]." ".$self->{resno}[$self->{resCA}[$i]].$self->{insertion_code}[$self->{resCA}[$i]]."\n";
	}
    }
    #}
    #else
    #{
	#my $len=scalar(@{$self->{resCA}});
	#for(my $i=0;$i<$len;$i++)
	#{
	#    push(@resno, $self->{resno}[$self->{resCA}[$i]].$self->{insertion_code}[$self->{resCA}[$i]]);
	#    print "då\n";
	#}
    #}
    #foreach my $elem (@resno)
    #{
#	print $elem."\n";
#    }
    return @resno;
}

=head2 print_pdb

    Title   : print_pdb
    Usage   : Prints out a pdb all atoms, where the SEQRES and
              CA residues match.
    Function: 
    Example : $pdb_obj->print_pdb(*FILEHANDLE, %lookup_hash)
    Returns : Void
    Args    : A filehandle  and a lookup_hash for the numbering of the residues (Optional);

=cut

sub print_pdb
{
    my ($self,$filehandle,%lookup_hash)=@_;
    my $use_lookup=scalar(keys(%lookup_hash));
    #my @ali = split(//,$self->{seqCA_ali});
    my @CAali = @{$self->{resCA}};
    		    
#print $self->{seqCA_ali};
    #print length($self->{seqCA_ali});
    #print "\n\n";
    #print $self->{seqCA}->seq();
    #print "\n\n";
    my $len = number_of_res($self);
    my $counter = 0;
    my $rows = 1;
    my $residues="";
    my $chain=$self->{chain};
    my @header=@{$self->{header}};
    if($use_lookup)
    {
	print $filehandle "REMARK The numbering of the residue changed to match the\nREMARK correct structure: ".scalar(localtime)."./BW\n";
    }
    #foreach my $line(@{$self->{header}})
    #print $#ali;
    #print "\n";
    #foreach my $line(@header)
    #{
#	print $filehandle $line."\n";
#    }
    
    for (my $i=0;$i<=$#CAali;$i++)
    {
	if($CAali[$i] ne "-")
	{
	    if($counter%13==0 && $counter != 0)
	    {
		
		printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
		#printf  ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
		$rows++;
		$residues="";
	    }
	    my $current_index=$self->{resCA}[$i];
	    #print $current_index."\n";
	    #print $i."\n";
	    $residues.=$self->{res}[$current_index]." ";
	    #print $residues."\n";
	    $counter++;
	}
    }
    if(length($residues) != 0)
    {
	printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
	#printf  ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
	$residues="";
	#$counter++;
    }
    
    if($use_lookup)
    {   	
	my $i=0;
	my $j=0;
	my $old_resname="undef";
	my $print_res=0;
	my $last_added_index=0;
	for($i=0;$i<scalar @{$self->{res}};$i++)
	{  
	    my $current_resname="$self->{resno}[$i]$self->{insertion_code}[$i]";
	    if($old_resname ne $current_resname)
	    {
		if($self->{resCA}[$j] eq "-")  # Just print the residue for which a CA exists
		{
		    $print_res=0  
		}
		else
		{
		    $print_res=1;
		}
		$j++;
	    }
	    
	    if($print_res)
	    { 
		printf $filehandle ("%-6s %4d  %-3s%1s%3s %1s%5s    %7.3lf %7.3lf %7.3lf\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$chain,$lookup_hash{$self->{resno}[$i].$self->{insertion_code}[$i]},$self->{x}[$i],$self->{y}[$i],$self->{z}[$i]);
		$last_added_index=$i;
	    }
	    $old_resname=$current_resname;
	}
	printf $filehandle ("%-6s %4d  %-3s%1s%3s %1s%5s\n","TER",$self->{atomno}[$last_added_index]+1,"" ,"" ,$self->{res}[$last_added_index],$chain,$lookup_hash{$self->{resno}[$last_added_index].$self->{insertion_code}[$last_added_index]});
	print $filehandle "END\n";
	close($filehandle);
    }
    else
    {
	my $i=0;
	my $j=0;
	my $old_resname="undef";
	my $print_res=0;
	my $last_added_index=0;
	for($i=0;$i<scalar @{$self->{res}};$i++)
	{  
	    #print "$self->{res}[$i] $self->{resno}[$i]$self->{insertion_code}[$i] $old_resname\n";
	    my $current_resname="$self->{resno}[$i]$self->{insertion_code}[$i]";
	    if($old_resname ne $current_resname)
	    {
		#print "In here $self->{resbegin}[$j]\n";

		if($self->{resCA}[$j] eq "-")
		{
		    $print_res=0;   # Changed to 1 from 0, which makes the whole if statement obsolete,
		}                   # I can't recall why I added it in the first place. In the case
		                    # CA is missing I think it shouldn't print the res.
		else
		{
		    $print_res=1;
		}
		$j++;
	    }
	    
	    if($print_res)
	    {
		printf $filehandle ("%-6s %4d  %-3s%1s%3s %1s%4d%1s   %8.3lf%8.3lf%8.3lf\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$chain,$self->{resno}[$i],$self->{insertion_code}[$i],$self->{x}[$i],$self->{y}[$i],$self->{z}[$i]);
		$last_added_index=$i;
	    }
		$old_resname=$current_resname;
	}
	printf $filehandle ("%-6s %4d  %-3s%1s%3s %1s%4d%1s\n","TER",$self->{atomno}[$last_added_index]+1,"" ,"" ,$self->{res}[$last_added_index],$chain,$self->{resno}[$last_added_index],$self->{insertion_code}[$last_added_index]);
	print $filehandle "END\n";
	close($filehandle);
    }
}

=head2 print_sliced_pdb

    Title   : print_sliced_pdb
    Usage   : Prints out a pdb all atoms, where the SEQRES and
              CA residues match in the interval specified.
    Function: 
    Example : $pdb_obj->print_pdb(*FILEHANDLE,$start,$end)
    Returns : Void
    Args    : A filehandle  and start and end points

=cut

sub print_sliced_pdb
{
    my ($self,$filehandle,$start,$end)=@_;
    #my @ali = split(//,$self->{seqCA_ali});
    #print $self->{seqCA_ali};
    #print length($self->{seqCA_ali});
    #print "\n\n";
    #print $self->{seqCA}->seq();
    #print "\n\n";
    my @CAali = @{$self->{resCA}}; #aligned to "first atom res"
    my $len = number_of_res($self);
    my $counter = 0;
    my $res_count=0;
    my $rows = 1;
    my $residues="";
    my $chain=$self->{chain};
    my @header=@{$self->{header}};
    if(($end-$start)<=$len || $start<=$end)
    {
	$len=$end-$start+1;
	print $filehandle "REMARK Only the residues from $start to $end /BW\n";
	#foreach my $line(@{$self->{header}})
	#print $#ali;
	#print "\n";
	foreach my $line(@header)
	{
	    print $filehandle $line."\n";
	}
	for (my $i=0;$i<=$#CAali;$i++)
	{
	    if($CAali[$i] ne "-" && $counter>=$start-1 && $counter<=$end-1)
	    {
		
		if($res_count%13==0 && $res_count != 0)
		{
		    
		    printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
		    #printf  ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
		    $rows++;
		    $residues="";
		}
		my $current_index=$self->{resCA}[$i];
		#print $current_index."\n";
		#print $i."\n";
		$residues.=$self->{res}[$current_index]." ";
		#print $residues."\n";
		$res_count++;
	    }
	    $counter++ if($CAali[$i] ne "-");
	}
	if(length($residues) != 0)
	{
	    printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
	    #printf  ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
	    $residues="";
	    #$counter++;
	}
	my $i=0;
	my $j=0;
	my $print_res=0;
	my $counter=-1;
	my $old_resname="";
	my $current_resname="";

	my $last_added_index=0;
	for($i=0;$i<scalar @{$self->{res}};$i++)
	{   
	    $current_resname="$self->{resno}[$i]$self->{insertion_code}[$i]";
	    if($current_resname ne $old_resname)
	    {
		if($self->{resCA}[$j] eq "-")
		{
		    $print_res=0
		}
		else
		{
		    $print_res=1;
		}
		$counter++;
	        $j++;
	    }

	    if(($counter>=$start-1 && $counter<=$end-1) && $print_res)
	    {
		printf $filehandle ("%-6s %4d  %-3s%1s%3s %1s%4d%1s   %8.3lf%8.3lf%8.3lf\n", $self->{field}[$i],$self->{atomno}[$i],$self->{atomtype}[$i],$self->{atom_alt_loc}[$i],$self->{res}[$i],$chain,$self->{resno}[$i],$self->{insertion_code}[$i],$self->{x}[$i],$self->{y}[$i],$self->{z}[$i]);
		$last_added_index=$i;
	    }
	    $old_resname=$current_resname;
	}
	printf $filehandle ("%-6s %4d  %-3s%1s%3s %1s%4d%1s\n","TER",$self->{atomno}[$last_added_index]+1,"" ,"" ,$self->{res}[$last_added_index],$chain,$self->{resno}[$last_added_index],$self->{insertion_code}[$last_added_index]);
	print $filehandle "END\n";
	close($filehandle);
    }
    else
    {
	warn "print_sliced_pdb: The starting point ($start) is after the end point ($end)\n" if($start>$end);
	warn "print_sliced_pdb: Specified interval $start-$end is larger than the length\n" if($end-$start>$len);
	return 0;
	    
    }
}





=head2 print_ca_pdb

    Title   : print_ca_pdb
    Usage   : Prints out a pdb with the CA atoms only, where the SEQRES and
              CA atom residues match.
    Function: 
    Example : $pdb_obj->print_pdb(*filehandle)
    Returns : Void
    Args    : A filehandle CHECK IT OUT NOT SURE IT IS FIXED

=cut


sub print_ca_pdb
{
    my ($self,$filehandle)=@_;
    my @ali = split(//,$self->{seqCA_ali});
    my $len = scalar @ali;
    my $counter = 0;
    my $rows = 1;
    my $residues="";
    my $chain=$self->{chain};
    
   #SEQRES
    if($self->{seqres} != 0)
    {
        for (my $i=0;$i<=$#ali;$i++)
	{
	    if($ali[$i] ne "-")
	    {
		if($counter%13==0 && $counter != 0)
		{
		    
		    printf ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
		    $rows++;
		    $residues="";
		}
		my $current_index=$self->{resCA}[$i];
		$residues.=$self->{res}[$current_index]." ";
		$counter++;
		#$residues.=$aa123{$res}." ";
		#$counter++;
	    }
	}
	if(length($residues) != 0)
	{
	    printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
	    $residues="";
	    #$counter++;
	}
    }
    for(my $i=0;$i<scalar @ali;$i++)
    {   
	if($ali[$i] ne "-")
	{
	    printf ("%-6s %4d  %-3s%1s%3s %1s%4d%1s   %8.3lf%8.3lf%8.3lf\n", $self->{field}[$self->{resCA}[$i]],$self->{atomno}[$self->{resCA}[$i]],$self->{atomtype}[$self->{resCA}[$i]],$self->{atom_alt_loc}[$self->{resCA}[$i]],$self->{res}[$self->{resCA}[$i]],$chain,$self->{resno}[$self->{resCA}[$i]],$self->{insertion_code}[$self->{resCA}[$i]],$self->{x}[$self->{resCA}[$i]],$self->{y}[$self->{resCA}[$i]],$self->{z}[$self->{resCA}[$i]]);
	}
    }
}






=head2 print_ca_model

    Title   : print_ca_model
    Usage   : Prints out a pdb with the CA atoms only, where the SEQRES and
              CA atom residues match the provided alignment
    Function: 
    Example : $pdb_obj->print_ca_model(*filehandle,$pdb_ali,$target_ali)
    Returns : Void
    Args    : A filehandle,alignment
    $pdb_ali sequence has to be taken from the correct pdbfile -> that sequence must match the exact CAres from the pdb

=cut


sub print_ca_model
{
    my ($self,$filehandle,$pdb_ali,$target_ali)=@_;
    
    my @index_to_CAres = @{$self->{resCA}};
    
    my @pdb_ali = split(//,$pdb_ali);
    my @target_ali = split(//,$target_ali);
    
    #print "$pdb_ali\n\n$target_ali\n\n";
    my $counter = 0;
    my $rows = 1;
    my $residues="";
    #my $chain=$self->{chain};
    my $chain=" ";
#count length,
    my $len = 0;
    for (my $i=0;$i<=$#pdb_ali;$i++)
    {
	$len++ if($pdb_ali[$i] ne '-' && $target_ali[$i] ne '-');
    }


#   #SEQRES
#    for (my $i=0;$i<=$#pdb_ali;$i++)
#    {
#	if($pdb_ali[$i] ne '-' && $target_ali[$i] ne '-')
#	{
#	    if($counter%13==0 && $counter != 0)
#	    {
#		
#		printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
#		$rows++;
#		$residues="";
#	    }
#	    $residues.=$aa123{$target_ali[$i]}." ";
#	    $counter++;
#	}
#    }
#    
#    if(length($residues) != 0)
#    {
#	printf $filehandle ("SEQRES  %2d %s %4d  %s\n", $rows,$chain,$len,$residues);
#	$residues="";
#	#$counter++;
#    }
    
    my $pdb_index=0;
    my $target_counter=0;
    for(my $target_index=0;$target_index<$#target_ali;$target_index++)
    {
	$target_counter++ if($target_ali[$target_index] ne '-');
	if($pdb_ali[$target_index] ne '-' && $target_ali[$target_index] ne '-')
	{
	    printf $filehandle ("%-6s %4d  %-3s%1s%3s %1s%4d%1s   %8.3lf%8.3lf%8.3lf\n", $self->{field}[$index_to_CAres[$pdb_index]],$self->{atomno}[$index_to_CAres[$pdb_index]],"CA",$self->{atom_alt_loc}[$index_to_CAres[$pdb_index]],$aa123{$target_ali[$target_index]},$chain,$target_counter," ",$self->{x}[$index_to_CAres[$pdb_index]],$self->{y}[$index_to_CAres[$pdb_index]],$self->{z}[$index_to_CAres[$pdb_index]]);
	    
	}
	$pdb_index++ if($pdb_ali[$target_index] ne '-');
    }
    printf $filehandle "TER\n";
}



=head2 fix_numbering

    Title   : fix_numbering
    Usage   : prints out a pdb file with the residue numbers (residue names) provided by the 
              correct pdb file

    Function: 
    Example : $pdb_obj->fix_numbering($correct_pdb_obj,*FILE);
    Returns : Void
    Args    : A correct pdb file, a filehandle

=cut

sub fix_numbering
{
    my ($model_obj,$correct_obj,$filehandle)=@_;
    my $model_seq=$model_obj->CAres()->seq();
    my $correct_seq=$correct_obj->CAres()->seq();
    my @resname_correct=$correct_obj->CAresname();
    my @resname_model=$model_obj->CAresname();
    
    my %lookup_hash=();
    #print scalar @resname_correct;
    #print "\n";
    #print length $correct_seq;
    #print "\n";
	
    my($model_seq2,$correct_seq2)=_align($model_seq,$correct_seq);
   # print "$model_seq2\n\n$correct_seq2\n\n";

    
    my @model_seq2=split('',$model_seq2);
    my @correct_seq2=split('',$correct_seq2);
    my $len=length($model_seq2);
    my $model_count=0;
    my $correct_count=0;
    my $missing_res_count=0;
    my $last_resnum=0;
    #my $special_counter=-999;
    my @temp=split(/[A-Z]/,$resname_model[$correct_count]);
    my $current_resnum=$temp[0];
    for(my $i=0;$i<$len;$i++)
    {
	if($correct_seq2[$i] ne "-" && $model_seq2[$i] ne "-")
	{
	    $lookup_hash{$resname_model[$model_count]}=$resname_correct[$correct_count];
	    my @temp=split(/[A-Z]/,$resname_correct[$correct_count]);
	   # print "$resname_correct[$correct_count] $resname_model[$model_count]\n";
	    $current_resnum=$temp[0];
	    $last_resnum=$current_resnum;
	    $model_count++;
	    $correct_count++;
	    #print "$resname_correct[$correct_count] $current_resnum\n";
	    #   print "$correct_seq2[$i] $model_seq2[$i] $resname_correct[$correct_count-1]\n";
	}
	elsif($correct_seq2[$i] ne "-")  # $model_seq2[$i] is not eq "-" because of first if.
	{
	  #  print $correct_seq2[$i]."\n";
	    
	    my @temp=split(/[A-Z]/,$resname_correct[$correct_count]);
	    $current_resnum=$temp[0];
	    $last_resnum=$current_resnum;
	    $correct_count++;
		
	}
	elsif($model_seq2[$i] ne "-")  # $correct_seq2[$i] is not eq "-" because of first if.
	{
	    $last_resnum++;
	    $missing_res_count=$last_resnum;
	    print "$missing_res_count"."X $resname_model[$model_count]\n";
	    $lookup_hash{$resname_model[$model_count]}="$missing_res_count"."X";
	    $model_count++;
	}
	
    }
    #foreach my $key (sort keys(%lookup_hash))
    #{
#	print "$key $lookup_hash{$key}\n";
#    }

    $model_obj->print_pdb($filehandle,%lookup_hash);

}

sub _align   # Takes two strings removes all dashes and returns the alignment. 
{
    my ($seq1,$seq2)=@_;
    $seq1=_remove_dashes($seq1);
    $seq2=_remove_dashes($seq2);
    #print $seq1."\n";
    #print $seq2."\n\n";
    my ($ali_return1,$ali_return2);
    my $factory=new Bio::Tools::pSW('-matrix' => '/home/ae/bjorn/modules/bioperl-live/examples/blosum62.bla','-gap' => 1,'-ext' => 0);
    my $seq_obj1=Bio::Seq->new( -seq => $seq1, -id => "seq1");
    my $seq_obj2=Bio::Seq->new( -seq => $seq2, -id => "seq2");

    my $aln = $factory->pairwise_alignment($seq_obj1,$seq_obj2);
    #$factory->align_and_show($seq_obj1,$seq_obj2,*STDOUT);
    ($ali_return1,$ali_return2)=_fix_alignment($aln,$seq_obj1,$seq_obj2);
    return ($ali_return1,$ali_return2);
}


sub distance
{
    my ($x1,$y1,$z1,$x2,$y2,$z2)=@_;
    #print "x1=$x1 y1=$y1 z1=$z1 x2=$x2 y2=$y2 z2=$z2\n";
    #print "\n";
    return sqrt(($x1-$x2)*($x1-$x2) + ($y1-$y2)*($y1-$y2) + ($z1-$z2)*($z1-$z2));
}

sub _remove_dashes
{
    my $temp=shift;
    my @list=split(/\-+/,$temp);
    my $new_str=join('',@list);
    return $new_str;

}



sub _fix_alignment
{
    my($aln,$seq1,$seq2)=@_;
    
    # Parse the alignment
    my ($seq,@start, @end,@ali_seqs);
    my $i=0;

    foreach $seq ($aln->each_seq())
    {
	$start[$i]=$seq->start();
	$end[$i]=$seq->end();
	$ali_seqs[$i]=$seq->seq();
	$i++;
    }
    #print $seq1->seq()."\n\n";
    #print $seq2->seq()."\n\n";
    #print $start[0]."\n".$start[1]."\n";
# Reformat alignment so that it contain all resides and -.
    
    my $ali_seq1="";
    my $ali_seq2="";
    
# Fix the begining
    if($start[0]!=1)
    {
	$ali_seq1.=$seq1->subseq(1,$start[0]-1);
	$ali_seq2.=_dashes($start[0]-1);
    }
    if($start[1]!=1)
    {
	$ali_seq1.=_dashes($start[1]-1);
	$ali_seq2.=$seq2->subseq(1,$start[1]-1);
    }
    # Add the alignment
    
    $ali_seq1.=$ali_seqs[0];
    $ali_seq2.=$ali_seqs[1];
    
# Fix alignment end;
    if($end[0]<$seq1->length())
    {
	my $len=$seq1->length();
	$ali_seq1.=$seq1->subseq($end[0]+1,$len);
	$ali_seq2.=_dashes($len-$end[0]);
    }
    
    if($end[1]<$seq2->length())
    {
	my $len=$seq2->length();
	$ali_seq1.=_dashes($len-$end[1]);
	$ali_seq2.=$seq2->subseq($end[1]+1,$len);
    }
    #print $ali_seq1."\n\n".$ali_seq2."\n\n";
    return ($ali_seq1,$ali_seq2);
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



