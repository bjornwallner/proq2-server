#!/usr/bin/perl -w

#input pdbfile.
#%param=@ARGV;
#$pdb=$param{-pdb};
use Cwd 'abs_path'; 
use File::Basename;
use File::Temp;
#`rm -fr /scratch/CASP9-TMP.*`;
my $TMP="/tmp/ProQ2.TMP.$$";
mkdir($TMP) if(!-e $TMP);
$caspdir="/afs/pdc.kth.se/home/b/bjornw/CASP9/";
$casp8dir="/afs/pdc.kth.se/home/b/bjornw/CASP9/casp8/";
$caspTMRdir="/afs/pdc.kth.se/home/b/bjornw/CASP9/TMR/";
$type=$ARGV[0]; #-fasta or -pdb
$full_path=abs_path($ARGV[1]);

my $script_path=dirname(abs_path($0));
#print $script_path."\n";                                                       
my $install_dir=abs_path("$script_path/../../../");
#print $install_dir."\n";        

#if(defined($ARGV[1]))
#{
#    if(not($ARGV[1]=~/overwrite/)) {
#	$source_folder=abs_path($ARGV[1]);
#    }
#}

#@temp=split(/\//,$full_path);
$path=dirname($full_path); #join('/',@temp[0..$#temp-1]);
#$path=abs_path($full_path); #join('/',@temp[0..$#temp-1];)
#print $path."\n";
#exit;
$casp_target=1;
$pdb=basename($full_path);




my $num="";
if($full_path=~/T0(\d\d\d)/ ||
    $full_path=~/TR(\d\d\d)/) {
    $target="T0$1";
    $num=$target;
    $num=~s/^T0//g;
    if($num<=514) {
	$caspdir=$casp8dir;
    }
    $targetdir=$caspdir;
    $casp_target=1;
} 
elsif($full_path=~/(TMR0\d)/) {
   
    $target=$1;
    $targetdir=$caspTMRdir;
    $caspdir=$caspTMRdir;
    $casp_target=1;
    print "TMR target, $target, look will for profiles here: $targetdir \n";
	
}
else
{
#    print STDERR "$full_path\n";
#    print STDERR "Not a CASP target you're in deep water now....I'll continue put on your lifevest! This will take time....\n";
    $targetdir=$path;
    my @temp=split(/\./,$pdb);
    $target=$temp[0];
    $casp_target=0;
}

print "casp_target=$casp_target\n";
chdir($path);
$refinement=0;
print "run all for $pdb in $path...\n";

if(defined($ARGV[2]))
{
    if($ARGV[2]=~/overwrite/) {

	print "Yes OVERWRITE!\n";
	`rm -f $pdb.chk`;
	`rm -f $pdb.mtx`;
	`rm -f $pdb.psi`;
	`rm -f $pdb.proqres2`;
	`rm -f $pdb.proqres2.gz`;
	`rm -f $pdb.rsa`;
	`rm -f $pdb.stride`;
	`rm -f $pdb.topcons`;
	`rm -f $pdb.topcons.fa`;
	`rm -f $pdb.ss2`;
	`rm -f $pdb.mpSA`;
	`rm -f $pdb.Z`;
	`rm -f $pdb.fasta`;
	`rm -f $pdb.zpred`;
	`rm -f $pdb.acc`;
	`rm -f $pdb.proq2`;
	
    }
    if($ARGV[2]=~/refinement/) {
	$refinement=1;
	
    }

}


my $subunits=`grep -c TER $pdb`;
chomp($subunits);
$subunits=1 if($subunits==0 || !$casp_target);

`mv -f $pdb $pdb.orig`;
`grep ^ATOM $pdb.orig | $install_dir/apps/ProQ2/bin/kill_chain.pl  > $pdb`;

#if(!$casp_target)
#fasta
$fasta="$pdb.fasta";
if($type eq "-fasta") {
    `cp $pdb $fasta`;
} else {
    $seq=`$install_dir/apps/ProQ2/bin/aa321CA.pl $pdb`;
    print $fasta."\n";
    if(!-e $fasta) {
	print "Creating fasta file..\n";
	
	open(OUT,">$fasta");
	print OUT ">$pdb\n$seq\n";
	close(OUT);
    }
}


$profile_file="$pdb.psi";
$profile_mtx_file="$pdb.mtx";
$ss2_file="$pdb.ss2";
#exit;
#
if($refinement) {

    $profile_file2="$caspdir/refinement/TR$num.pdb.psi";
    $profile_mtx_file2="$caspdir/refinement/TR$num.pdb.mtx";
    $ss2_file2="$caspdir/refinement/TR$num.pdb.ss2";
    symlink($profile_file2,$profile_file);
    symlink($profile_mtx_file2,$profile_mtx_file);
    symlink($ss2_file2,$ss2_file);
}

if(!-e $profile_file || !-e $profile_mtx_file || !-e $ss2_file) 
{
    print "Creating profiles and predicting ss\n";
    my $tdir=`pwd`;
    print $tdir."\n";
   # print "/Users/bjornw/Research/bin/create_profile.sh $fasta\n";
   # `/Users/bjornw/Research/bin/create_profile.sh $fasta`;

    
    $target_profile=$caspdir."/$target/$target.psi";
    $target_profile_mtx=$caspdir."/$target/$target.mtx";
    $target_ss2=$caspdir."/$target/$target.ss2";

    if(!$casp_target) {
	$target_profile="$targetdir/$target.psi";
	$target_profile_mtx="$targetdir/$target.mtx";
	$target_ss2="$targetdir/$target.ss2";

    }
    
    #print "/Library/WebServer/ProQM/bin/create_profile.sh $fasta\n";
    #print " 
    if($casp_target && (!-e $target_profile ||
			!-e $target_profile_mtx || 
			!-e $target_ss2))
    {
	print STDERR "Need to have prediction for full length sequence! run init_targets_elofsson.pl!\n";
	print STDERR "$target_profile\n$target_profile_mtx\n$target_ss2\n";
	exit
    }
    if(!-e $target_profile ||
       !-e $target_profile_mtx || 
       !-e $target_ss2)
    {
	print STDERR "Running PSI-BLAST (this will also create some input to accpro...)\n";
	print "$install_dir/bin/create_profile.sh $fasta\n";
	`$install_dir/bin/create_profile.sh $fasta`;
    } else { 
   
	print "$install_dir/apps/ProQ2/bin/profile_subset.pl $target_profile $fasta $profile_file\n";
	#exit;
	system("$install_dir/apps/ProQ2/bin/profile_subset.pl $target_profile $fasta $profile_file $subunits");
	system("$install_dir/apps/ProQ2/bin/profile_subset.pl $target_profile_mtx $fasta $profile_mtx_file $subunits");
	system("$install_dir/apps/ProQ2/bin/ss2_subset.pl $target_ss2 $fasta $ss2_file $subunits");
    }

}
#exit;
if(!-e "$pdb.stride") {
    print "stride\n";
    `$install_dir/apps/stride/linux/bin/stride $pdb.orig >$pdb.stride`;
}


$naccessfile="$pdb.rsa";
if(!-e $naccessfile) {
    print "Naccess\n";
    $tmpfile=$$;
    chdir($TMP);
    `ln -s $path/$pdb $tmpfile.pdb`;
    #print "ln -s $path/$pdb $tmpfile.pdb\n";
    #print "$tmpfile\n";
    #exit;
    `$install_dir/apps/naccess/naccess $tmpfile.pdb`;
    print "$install_dir/apps/naccess/naccess $tmpfile.pdb\n";
#print "$tmpfile.rsa\n";
    if(1==0) {
	open(RSA,"$tmpfile.rsa");
	open(OUT,">$path/$naccessfile");
	my $old_chain="undef";
	while(<RSA>) {
	    if(/^RES/) {
		if($old_chain eq "undef") {
		    $old_chain=substr($_,8,1);
		    
		}
		$chain=substr($_,8,1);
		if($chain ne $old_chain) {
		    last;
		} else
		{
		    print OUT;
		}
	    }
	    
	}
	close(RSA);
    }
    `cp $tmpfile.rsa $path/$naccessfile`;
    
   # exit;
    `rm -fr $tmpfile.*`;
    chdir($path);
}

$accfile="$pdb.acc";

if($refinement) {
    $accfile2="$caspdir/refinement/TR$num.pdb.acc";
    symlink($accfile2,$accfile);
}

if(!-e $accfile) {
    print "Accpro\n";
    $target_acc=$targetdir."/$target/$target.acc";
    if(!$casp_target) {
	$target_acc=$targetdir."/$target.acc";
    }
    if(!-e $target_acc) {

	#system("$install_dir/apps/sspro4/bin/predict_acc.sh $fasta $accfile");
	`$install_dir/apps/sspro4/bin/predict_acc.sh $fasta $accfile`;
    } else {
	print "subset with $target_acc\n";
	system("$install_dir/apps/ProQ2/bin/acc_subset.pl $target_acc $fasta $accfile $subunits");
    }
    
}


$proqres_flags="-rc 6 -ac 4 -w 23 -t 6";
$proqoutfile="$pdb.proqres2";
$proqoutfile_gz="$pdb.proqres2.gz";
#print "/afs/pdc.kth.se/home/b/bjornw/source/c/pdb/ProQres/PROQRES/bin/ProQres.Darwin -pdb $pdb -surf $naccessfile -stride $pdb.stride -psipred $pdb.ss2 -prof $pdb.psi $proqres_flags -output_input\n";

if(-e $proqoutfile_gz) {

    if(-s $proqoutfile_gz < 100) {

	`rm $proqoutfile_gz`;
    }
}

if(!-e $proqoutfile_gz) 
{
    print "ProQres...\n";
    
    if(-e "/Users/bjorn/") {
	`export PROQRESDIR=$install_dir/apps/ProQres/weights; $install_dir/apps/ProQres/bin/ProQres -pdb $pdb -surf $naccessfile -stride $pdb.stride -psipred $pdb.ss2 -prof $pdb.psi $proqres_flags -output_input > $TMP/$proqoutfile.$$`;
    } else {

	print "export PROQRESDIR=$install_dir/apps/ProQres/weights; $install_dir/apps/ProQres/bin/ProQres -pdb $pdb -surf $naccessfile -stride $pdb.stride -psipred $pdb.ss2 -prof $pdb.psi $proqres_flags -output_input\n";
	`ulimit -s unlimited;export PROQRESDIR=$install_dir/apps/ProQres/weights; $install_dir/apps/ProQres/bin/ProQres64 -pdb $pdb -surf $naccessfile -stride $pdb.stride -psipred $pdb.ss2 -prof $pdb.psi $proqres_flags -output_input > $TMP/$proqoutfile.$$`;
    }
    `gzip $TMP/$proqoutfile.$$`;
    `mv $TMP/$proqoutfile.$$.gz $proqoutfile.gz`;
}
 
`rm -fr $TMP`;
