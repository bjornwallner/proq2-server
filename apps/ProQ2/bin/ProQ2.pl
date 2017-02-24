#!/usr/bin/perl -w
use strict;
use File::Basename;
use Cwd 'abs_path'; 
use Digest::MD5 qw(md5 md5_hex md5_base64);
#exit
# This is a fire and forget implementation of ProQ2
# that will run all billions of external programs
# and produce a fantastic result!!!!

# all files will be produced in the directory where the pdb lies
my $args=join(" ",@ARGV);
print "CMD: ProQ2.pl $args\n";
my $script_path=dirname(abs_path($0));
#print $script_path."\n";
my $install_dir=abs_path("$script_path/../../../");
#print $install_dir."\n";
#exit;
my $verbose=1;
my $cleanup=0;
my $SVM_MODEL_DIR="$install_dir/apps/ProQ2/data/";
my $SVM_CLASSIFY="$install_dir/apps/svm_light/svm_classify";
my $pdb=$ARGV[0];
my $email_file=$ARGV[1];
my $seqres_file=$ARGV[2];
my $bfactorPDB=0;
if(defined($ARGV[3])) {
    $bfactorPDB=1;
    print "Will output bfactorPDB!\n";
}

my $sendemail=1;
 
my $email=`cat $email_file | head`;
chomp($email);

$sendemail=0 if($email_file=~/noemail/ || 
		$email eq "noemail" ||
		$email=~/fake/);

my $folder=dirname(abs_path($pdb)); #join("/",@temp[0..$#temp-1]);
my $seq_cache="$install_dir/server/output/cache/";
`mkdir $seq_cache` if(!-e $seq_cache);
my $pdb_base=basename($pdb); #$temp[$#temp];
my @temp=split(/\//,$folder);
my $jobid=$temp[$#temp];
my $subfolder=$temp[$#temp-1];
my $job_folder="$subfolder/$jobid";
print $folder."\n";


my $outfile="$folder/$pdb_base.proq2";
my  $target_length=`grep -c ' CA ' $pdb`;
chomp($target_length);
my $subunits=`grep -c TER $pdb`;
$subunits=1 if($subunits==0);

#before launching cache the seqres sequence or pdb.....
#get the seqres

my $seqres="";
open(FILE,$seqres_file);
while(<FILE>) {
    chomp;
    if(not(/^/)) {
	$seqres.=$_;
    }
}
$seqres=~s/\s+//g;
#$seqres=~s/\n//g;
my $pdbseq=`$install_dir/apps/ProQ2/bin/aa321CA.pl $pdb`;
open(OUT,">$pdb.fasta");
print OUT ">$pdb\n$pdbseq\n";
close(OUT);
my $md5pdbseq = md5_hex($pdbseq);
my $pdbseq_folder="$seq_cache/$md5pdbseq";
my $md5seqres = md5_hex($seqres);
my $seqres_folder="$seq_cache/$md5seqres";
my $cache_after=0;
if(length($seqres)>10) {
    if(!-e $seqres_folder) {
	`mkdir $seqres_folder`;
	open(OUT,">$seqres_folder/seqres");
	print OUT ">seqres\n$seqres\n";
	close(OUT);
	`$install_dir/apps/ProQ2/bin/run_all_external.pl -fasta $seqres_folder/seqres`;
	print "Caching seqres...\n";
    } else {
	print "seqres already in cache $seqres_folder!\n";
    }
    #collect subset data
#    my $PERL5LIB="export PERL5LIB=/local/www/services/ProQ2/apps/ProQ2/bin";
 #   `$PERL5LIB export PERL5LIB=/local/www/services/ProQ2/apps/ProQ2/bin;`;
    print "getting data from cache...\n";
    `$install_dir/apps/ProQ2/bin/acc_subset.pl $seqres_folder/seqres.acc $pdb.fasta $pdb.acc`;
    `$install_dir/apps/ProQ2/bin/profile_subset.pl $seqres_folder/seqres.psi $pdb.fasta $pdb.psi`;
    `$install_dir/apps/ProQ2/bin/profile_subset.pl $seqres_folder/seqres.mtx $pdb.fasta $pdb.mtx`;
    `$install_dir/apps/ProQ2/bin/ss2_subset.pl $seqres_folder/seqres.ss2 $pdb.fasta $pdb.ss2`;
} elsif(-e $pdbseq_folder) {
    print "pdbseq already in cache $pdbseq_folder\n";
    `gunzip -f $pdbseq_folder/*.gz`;
    `cp $pdbseq_folder/seqres.acc $pdb.acc`;
    `cp $pdbseq_folder/seqres.psi $pdb.psi`;
    `cp $pdbseq_folder/seqres.mtx $pdb.mtx`;
    `cp $pdbseq_folder/seqres.ss2 $pdb.ss2`;
    `gzip -f $pdbseq_folder/*`;
} else {
        $cache_after=1;
}

if($verbose) {
    system("$install_dir/apps/ProQ2/bin/run_all_external.pl -pdb $pdb");
} else {
    `$install_dir/apps/ProQ2/bin/run_all_external.pl -pdb $pdb`;# -overwrite`;
}

if($cache_after) {
    print "Caching pdbseq...\n";
    `mkdir $pdbseq_folder`;
    `cp $pdb.acc $pdbseq_folder/seqres.acc`;
    `cp $pdb.psi $pdbseq_folder/seqres.psi`;
    `cp $pdb.mtx $pdbseq_folder/seqres.mtx`;
    `cp $pdb.ss2 $pdbseq_folder/seqres.ss2`;
    `gzip -f $pdbseq_folder/*`;
}


#print "/local/www/services/ProQ2/apps/ProQ2/bin/generate_svm_input.pl -i $pdb -o $folder/ -classify 1 -noprefix 1 -Sscore 3 -rc 6 -ac 4 -atom 1 -res 1 -surf50 1 -surf100 1 -surf25 1 -surf75 1 -pw 1 -pwin 23 -stride 5 -ss 1 -ss_sc 21 -entropy 3 -profile 0 -termini 5 -rsa_sc 21 -rsa 13 -grsa_sc 1 -gss_sc 1\n";

print "$install_dir/apps/ProQ2/bin/generate_svm_input.pl -i $pdb -o $folder/ -classify 1 -noprefix 1 -Sscore 3 -rc 6 -ac 4 -atom 1 -res 1 -surf50 1 -surf100 1 -surf25 1 -surf75 1 -pw 1 -pwin 23 -stride 5 -ss 1 -ss_sc 21 -entropy 3 -profile 0 -termini 5 -rsa_sc 21 -rsa 13 -grsa_sc 1 -gss_sc 1\n";
`export PERL5LIB=$install_dir/apps/ProQ2/bin;$install_dir/apps/ProQ2/bin/generate_svm_input.pl -i $pdb -o $folder/ -classify 1 -noprefix 1 -Sscore 3 -rc 6 -ac 4 -atom 1 -res 1 -surf50 1 -surf100 1 -surf25 1 -surf75 1 -pw 1 -pwin 23 -stride 5 -ss 1 -ss_sc 21 -entropy 3 -profile 0 -termini 5 -rsa_sc 21 -rsa 13 -grsa_sc 1 -gss_sc 1`;
#system("/local/www/services/ProQ2/apps/ProQ2/bin/generate_svm_input.pl -i $pdb -o $folder/ -classify 1 -noprefix 1 -Sscore 3 -rc 6 -ac 4 -atom 1 -res 1 -surf50 1 -surf100 1 -surf25 1 -surf75 1 -pw 1 -pwin 23 -stride 5 -ss 1 -ss_sc 21 -entropy 3 -profile 0 -termini 5 -rsa_sc 21 -rsa 13 -grsa_sc 1 -gss_sc 1");


my $svm_para_file="$folder/svmdata-Sscore3-ac4-atom1-entropy3-grsa_sc1-gss_sc1-profile0-prsa0-pw1-pwin23-rc6-res1-restypes6-rsa13-rsa_sc21-ss1-ss_sc21-stride5-surf1001-surf251-surf501-surf751-termini5-topology0-z0-z_sc0-zpred0/$pdb_base.svm";
#svmdata-Sscore3-ac4-atom1-entropy11-profile3-prsa1-pw1-pwin21-rc6-res1-restypes6-rsa1-rsa_sc21-ss1-ss_sc21-stride11-surf1001-surf251-surf501-surf751-termini23-topology11-z1-z_sc0-zpred1/$pdb_base.svm";
if(!-e $svm_para_file) {

    print STDERR "Something went wrong skipping... no worries the model probably suck!\n";
    exit;
}
my $svm_pred_base=$svm_para_file.".pred"; 
#if(1==0) {
for(my $i=1;$i<=5;$i++) {

#    print "Classifying using model.$i ...\n";
    
#    print "$SVM_CLASSIFY $svm_para_file $SVM_MODEL_DIR/model.linear.$i $svm_pred_base.$i\n";
    `$SVM_CLASSIFY $svm_para_file $SVM_MODEL_DIR/model.linear.$i $svm_pred_base.$i`;
    
}
#}

my %pred=();
for(my $i=1;$i<=5;$i++) {

    #print "$svm_pred_base.$i\n";
    open(FILE,"$svm_pred_base.$i");
    while(<FILE>) {
	chomp;
	push(@{$pred{$i}},$_);
    }
    close(FILE);
}




my @mean=();
my @std=();
my @dist_dev=();
my $d0=3;
for(my $j=0;$j<scalar(@{$pred{1}});$j++) {
    my @vec=();
    for(my $i=1;$i<=5;$i++) {
	push(@vec,$pred{$i}[$j]);

    }
    my $m=mean(@vec);
    my $d=15;
    $m=1 if($m>1);
    $m=0 if($m<0);
    if($m>0.03846) {
	$d=sqrt(1/$m-1)*$d0;
    }
    push(@mean,$m);
    push(@dist_dev,$d);
    push(@std,std(@vec));
    

}
my $global_quality=sum(@mean);
my $global_quality_norm=$global_quality/$target_length;


#changeB($pdb.orig,@dist_dev)
open(OUT,">$pdb.orig.B");
open(IN,"$pdb.orig");
my $index=-1;
my $old_resnum="undef";
while(<IN>) {
    if(/^ATOM/) {
#       my $atomno=substr($_, 7, 4);
#       my $atomtypeq=substr($_, 13, 3);
        my $resnum=substr($_,22,5);
#        $resnum=~s/\s+//g;
        #print "$resnum $old_resnum $atomtype\n";
        if($old_resnum ne $resnum)
        {
            $index++;
#           print "POS $ali_resnum_model[$pos] $_";
        }
        $old_resnum=$resnum;
	substr($_,60,6)=sprintf("%6.2f",$dist_dev[$index]);
    }
    print OUT;


}
close(OUT);
close(IN);

my @CAs=`grep CA $pdb.orig`;
my @id=();
foreach my $ca(@CAs) {
    my $chain=substr($ca,21,1);
    my $resno=substr($ca, 22, 4);
    my $id="$resno:$chain";
    push(@id,$id);
}


if($subunits>1) {

    print "sub $subunits\n";
    my @mean2=();
    for(my $i=0;$i<scalar(@mean)/$subunits;$i++) {
	push(@mean2,$mean[$i]);
    }
    $global_quality=sum(@mean2);
    $global_quality_norm=$global_quality/$target_length;
    
}


#my $email_body="Predicted quality [0,1] , 0-worst, 1-best\n\n";
#$email_body.="Global Quality: $global_quality\n";
#$email_body.="Local Quality:\n";
open(OUT,">$outfile");
print OUT "# first column pos, second column average of five last. third column std of five last\n";
for(my $j=0;$j<scalar(@{$pred{1}});$j++) {
    printf OUT ("%-5d %10.7f %10.7f ",$j+1,$mean[$j],$std[$j]);
    #$email_body.=sprintf("%-5d %-10.7f\n",$j+1,$mean[$j]);
    for(my $i=1;$i<=5;$i++) {
	printf OUT ("%10.7f ",$pred{$i}[$j]);

    }
    printf OUT ("\n");
}
printf OUT ("Global quality: %10.7f sum= %10.7f L= %4d subunits= %2d (Quality only calculated for the first subunit\n",$global_quality_norm,$global_quality,$target_length,$subunits);
close(OUT);
print "Output printed to $outfile\n";


my $email_bodyB="REMARK Link to predicted distance deviation plot: $install_dir/output/$job_folder/local.D.png\n";
$email_bodyB.="REMARK Raw files: http://duffman.it.liu.se/ProQ2/output/$job_folder/\n";
$email_bodyB.="REMARK Predicted local quality: [0,1] , 0-worst, 1-best, dist=3*sqrt(1/S-1)\n";
$email_bodyB.="REMARK Predicted global quality: sum of local S\n";
$email_bodyB.="REMARK Global Quality: ";

$email_bodyB.=sprintf("%-8.3f\n",$global_quality);
$email_bodyB.="REMARK Local predicted distance deviations are in the B factor column\n";
$email_bodyB.=`cat $pdb.orig.B`;
my $email_body="";
$email_body.="Link to predicted distance deviation plot: $install_dir/output/$job_folder/local.D.png\n";
#$email_body.="Link to local quality plot: http://proq2.theophys.kth.se/output/$jobid/local.png\n";
$email_body.="Raw files: http://duffman.it.liu.se/ProQ2/output/$job_folder/\n";
$email_body.="Predicted local quality: [0,1] , 0-worst, 1-best, dist=3*sqrt(1/S-1)\n";
$email_body.="Predicted global quality: sum of local S\n";
$email_body.="\n";

$email_body.="Global Quality: ";
$email_body.=sprintf("%-8.3f\n",$global_quality);
$email_body.="Local Quality:\n";
$email_body.=sprintf("%-5s %-8s %-8s %-8s\n","Pos","S","dist","resnum:chain");
open(OUT,">$svm_pred_base.all");
print OUT "# first column pos, second column average of five last. third distance, fourth column raw average, fifth column std of five last, fourth column resno:chain\n";
#my $d0=3;
for(my $j=0;$j<scalar(@{$pred{1}});$j++) {

    my $m=$mean[$j];
    $m=1 if($m>1);
    $m=0 if($m<0);
    my $d=15;
    if($m>0.03846) {
	$d=sqrt(1/$m-1)*$d0;
    }
    printf OUT ("%-5d %-10.7f %-10.7f %-10.7f %-10.7f %-8s ",$j+1,$m,$d,$mean[$j],$std[$j],$id[$j]);
    $email_body.=sprintf("%-5d %-8.3f %-8.3f %-8s\n",$j+1,$m,$d,$id[$j]);
    for(my $i=1;$i<=5;$i++) {
	printf OUT ("%-10.7f ",$pred{$i}[$j]);

    }
    printf OUT ("\n");
}

close(OUT);
print "Output printed to $svm_pred_base.all\n";

open(OUT,">$svm_pred_base.all.gnu");
print OUT "set terminal png font times\n";
print OUT "set output '$svm_pred_base.png'\n";
#print OUT "set multiplot layout 1,2\n";
print OUT "set xlabel 'sequence position'\n";
print OUT "set ylabel 'predicted quality, 1-best 0-worst'\n";
print OUT "plot '$svm_pred_base.all' u 1:2 w l lw 1 title 'raw','' u 1:2 w l lw 2 smooth bezier title 'smoothed'\n";

print OUT "set ylabel 'predicted distance deviation (angstroms)'\n";
print OUT "set output '$svm_pred_base.D.png'\n";
print OUT "plot '$svm_pred_base.all' u 1:3 w l lw 1 title 'raw','' u 1:3 w l lw 2 smooth bezier title 'smoothed'\n";


close(OUT);
system("gnuplot $svm_pred_base.all.gnu");
`cp $svm_pred_base.png $folder/local.png`;
`cp $svm_pred_base.D.png $folder/local.D.png`;

if($sendemail) {
#send email with result.
    
    $email=`cat $email_file | head`;
    my $target_name=`cat $folder/name|head`;
    
    chomp($email);
    chomp($target_name);
    if(length($target_name)>0) {
	$target_name="for $target_name";
    }
    my $extra="";
#    if($email=~/proteinmodel/) {
#	$extra=",bjornw\@ifm.liu.se";
#   }
    open(EMAIL,">$folder/results.email");
    print EMAIL "From: bjornw\@ifm.liu.se\n";
    print EMAIL "Subject: Local Quality Prediction ProQ2-server Results $target_name\n";
    print EMAIL "To: $email$extra\n";
    if($bfactorPDB) {
	print "Bfactor\n";
	print EMAIL $email_bodyB;
    } else {
	print "no Bfactor\n";
       	print EMAIL $email_body;
    }
    close(EMAIL);
    
    

    `cat $folder/results.email | /usr/sbin/sendmail -t`;
} else {

    `cp  $svm_pred_base.all $pdb.ProQ2`;
}


`cd $folder;tar -cf sequence_spec_input.tar $pdb_base.acc $pdb_base.fasta $pdb_base.mtx $pdb_base.psi $pdb_base.ss2`;

if($cleanup==1) {

    print "Cleanup...\n";
    unlink("$folder/$pdb_base.acc");
    unlink("$folder/$pdb_base.fasta");
    unlink("$folder/$pdb_base.mtx");
    unlink("$folder/$pdb_base.psi");
    unlink("$folder/$pdb_base.rsa");
    unlink("$folder/$pdb_base.ss2");
    unlink("$folder/$pdb_base.stride");
    unlink("$folder/$pdb_base.proqres2.gz");
    unlink($svm_para_file);
}
#
#open(OUT,">$folder/pred.all.gnu");
#print OUT "set terminal png font times\n";
#print OUT "set output '$folder/pred.png'\n";
#print OUT "set xlabel 'sequence position'\n";
#print OUT "set ylabel 'predicted quality, 1-best 0-worst'\n";
#print OUT "plot '$folder/pred.all' u 1:2 w l lw 1 title 'raw','' u 1:2 w l lw 2 smooth bezier title 'smoothed'\n";
#close(OUT);
#system("/opt/local/bin/gnuplot $folder/pred.all.gnu");
#
##send email with result.
#
#my $email=`cat $email_file | head`;
#my $target_name=`cat $folder/name|head`;
#
#chomp($email);
#chomp($target_name);
#if(length($target_name)>0) {
#    $target_name="for $target_name";
#}
#open(EMAIL,">$folder/results.email");
#print EMAIL "From: ProQM\@cbr.su.se\n";
#print EMAIL "Subject: Local Quality Prediction ProQM-server Results $target_name\n";
#print EMAIL "To: $email\n";
#print EMAIL $email_body;
#close(EMAIL);
#
#
#
#`cat $folder/results.email | /usr/sbin/sendmail -t`;


sub std
{
    my @data=@_;
    my $mean=mean(@data);
    my $n=scalar(@data);
    if($n!=1)
    {
	my $sum=0;
	foreach my $term(@data)
	{
	    $sum+=($term-$mean)*($term-$mean);
	}
	return sqrt(1/($n-1)*$sum);
    }
    else
    {
	return 0;
    }
}

sub mean
{
    my @data=@_;
    my $sum=0;
    my $number_of_elements=scalar @data;
    foreach my $term(@data)
    {
	$sum+=$term;
    }
    if($number_of_elements==0)
    {
	return 0;
    }
    else
    {
	return $sum/$number_of_elements;
    }
}


sub sum
{
    my @data=@_;
    my $sum=0;
    my $number_of_elements=scalar @data;
    foreach my $term(@data)
    {
	$sum+=$term;
    }
    if($number_of_elements==0)
    {
	return 0;
    }
    else
    {
	return $sum;
    }
}
