#!/usr/bin/perl -w

#input pdbfile.
#%param=@ARGV;
#$pdb=$param{-pdb};
use File::Basename;
use File::Temp;
$caspdir="/afs/pdc.kth.se/home/b/bjornw/CASP9/";
$casp8dir="/afs/pdc.kth.se/home/b/bjornw/CASP9/casp8/";
$full_path=$ARGV[0];
#@temp=split(/\//,$full_path);
$path=dirname($full_path); #join('/',@temp[0..$#temp-1]);
$pdb=basename($full_path);
if($full_path=~/(T0\d\d\d)/) {
    $target=$1;
    my $num=$target;
    $num=~s/^T0//g;
    if($num<=514) {
	$caspdir=$casp8dir;
    }
} 
else
{
    print STDERR "Need to have casp target T0XXX in the path somewhere...\n";
    exit;
}
my $host=`hostname`;
chomp($host);

chdir($path);
print "run all for $pdb in $path...\n";
if(defined($ARGV[1]))
{
    if($ARGV[1] =~/overwrite/) {

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

}


#fasta
$fasta="$pdb.fasta";
$seq=`/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/ProQ2/bin/aa321CA.pl $pdb`;
if(!-e $fasta) {
    print "Creating fasta file..\n";
    
    open(OUT,">$fasta");
    print OUT ">$pdb\n$seq\n";
    close(OUT);
}


$profile_file="$pdb.psi";
$profile_mtx_file="$pdb.mtx";
$ss2_file="$pdb.ss2";
#exit;
#

if(!-e $profile_file || !-e $profile_mtx_file || !-e $ss2_file) 
{
    print "Creating profiles and predicting ss\n";

   # print "/Users/bjornw/Research/bin/create_profile.sh $fasta\n";
   # `/Users/bjornw/Research/bin/create_profile.sh $fasta`;
    $target_profile=$caspdir."/$target/$target.psi";
    $target_profile_mtx=$caspdir."/$target/$target.mtx";
    $target_ss2=$caspdir."/$target/$target.ss2";
    
    
    #print "/Library/WebServer/ProQM/bin/create_profile.sh $fasta\n";
    #print " 
    if(!-e $target_profile ||
       !-e $target_profile_mtx || 
       !-e $target_ss2)
    {
	print STDERR "Need to have prediction for full length sequence! run init_targets_elofsson.pl!\n";
	exit
    }
    system("/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/ProQ2/bin/profile_subset.pl $target_profile $fasta $profile_file");
    system("/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/ProQ2/bin/profile_subset.pl $target_profile_mtx $fasta $profile_mtx_file");
    system("/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/ProQ2/bin/ss2_subset.pl $target_ss2 $fasta $ss2_file");

}
#exit;
if(!-e "$pdb.stride") {
    print "stride\n";
    `/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/stride/linux/bin/stride $pdb >$pdb.stride`;
}


$naccessfile="$pdb.rsa";
if(!-e $naccessfile) {
    print "Naccess\n";
    $tmpfile=$$;
    `ln -s $pdb $tmpfile.pdb`;
    if(not($host=~/elofsson/)) {
	`/afs/pdc.kth.se/home/b/bjornw/modules/naccess $tmpfile.pdb`;
    } else {
	`/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/naccess/naccess $tmpfile.pdb`;
    }
#print "$tmpfile.rsa\n";
    `cp $tmpfile.rsa $naccessfile`;
   # print "$tmpfile\n";
    `rm -fr $tmpfile.*`;
}

$accfile="$pdb.acc";
if(!-e $accfile) {
    print "Accpro\n";
    $target_acc=$caspdir."/$target/$target.acc";
    `/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/ProQ2/bin/acc_subset.pl $target_acc $fasta $accfile`;
    
}


$proqres_flags="-rc 6 -ac 4 -w 23 -t 6";
$proqoutfile="$pdb.proqres2";
#print "/afs/pdc.kth.se/home/b/bjornw/source/c/pdb/ProQres/PROQRES/bin/ProQres.Darwin -pdb $pdb -surf $naccessfile -stride $pdb.stride -psipred $pdb.ss2 -prof $pdb.psi $proqres_flags -output_input\n";
if(!-e "$proqoutfile.gz") 
{
    print "ProQres...\n";
  #  print "export PROQRESDIR=/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/ProQres/weights; /afs/pdc.kth.se/home/b/bjornw/CASP9/apps/ProQres/bin/ProQres64 -pdb $pdb -surf $naccessfile -stride $pdb.stride -psipred $pdb.ss2 -prof $pdb.psi $proqres_flags -output_input\n";
    `export PROQRESDIR=/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/ProQres/weights; /afs/pdc.kth.se/home/b/bjornw/CASP9/apps/ProQres/bin/ProQres64 -pdb $pdb -surf $naccessfile -stride $pdb.stride -psipred $pdb.ss2 -prof $pdb.psi $proqres_flags -output_input > $proqoutfile`;
    `gzip $proqoutfile`;
}
 
