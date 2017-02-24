#!/usr/bin/perl -w
use strict;
use File::Basename;
# This is a fire and forget implementation of ProQM 
# that will run all billions of external programs
# and produce a fantastic result!!!!

# all files will be produced in the directory where the pdb lies
# ideally this should be a uniq directory that will live for 30 days....
my $SVM_MODEL_DIR="/Library/WebServer/ProQM/SVM-MODELS/svmdata-Sscore3-ac4-atom1-entropy11-profile3-prsa1-pw1-pwin21-rc6-res1-restypes6-rsa1-rsa_sc21-ss1-ss_sc21-stride11-surf1001-surf251-surf501-surf751-termini23-topology11-z1-z_sc0-zpred1/";
my $SVM_CLASSIFY="/Library/WebServer/ProQM/bin/svm_classify";
my $pdb=$ARGV[0];
my $email_file=$ARGV[1];
my $sendemail=1;

$sendemail=0 if($email_file=~/noemail/);
my @temp=split(/\//,$pdb);
my $folder=dirname($pdb); #join("/",@temp[0..$#temp-1]);
my $pdb_base=basename($pdb); #$temp[$#temp];
my $jobid=$temp[$#temp-1];
print $folder."\n";


print "/Library/WebServer/ProQM/bin/run_all_external.pl $pdb\n";
#system("/Library/WebServer/ProQM/bin/run_all_external.pl $pdb");
`/Library/WebServer/ProQM/bin/run_all_external.pl $pdb`;
print "/Library/WebServer/ProQM/bin/generate_svm_input.pl -i $pdb -o $folder/ -classify 1 -Sscore 3 -ac 4 -atom 1 -entropy 11 -profile 3 -prsa 1 -pw 1 -pwin 21 -rc 6 -res 1 -restypes 6 -rsa 1 -rsa_sc 21 -ss 1 -ss_sc 21 -stride 11 -surf100 1 -surf25 1 -surf50 1 -surf75 1 -termini 23 -topology 11 -z 1 -z_sc 0 -zpred 1\n";
system("/Library/WebServer/ProQM/bin/generate_svm_input.pl -i $pdb -o $folder/ -classify 1 -noprefix 1 -Sscore 3 -ac 4 -atom 1 -entropy 11 -profile 3 -prsa 1 -pw 1 -pwin 21 -rc 6 -res 1 -restypes 6 -rsa 1 -rsa_sc 21 -ss 1 -ss_sc 21 -stride 11 -surf100 1 -surf25 1 -surf50 1 -surf75 1 -termini 23 -topology 11 -z 1 -z_sc 0 -zpred 1");


my $svm_para_file="$folder/svmdata-Sscore3-ac4-atom1-entropy11-profile3-prsa1-pw1-pwin21-rc6-res1-restypes6-rsa1-rsa_sc21-ss1-ss_sc21-stride11-surf1001-surf251-surf501-surf751-termini23-topology11-z1-z_sc0-zpred1/$pdb_base.svm";
my $svm_pred_base="$folder/pred.$$";
#if(1==0) {
for(my $i=1;$i<=5;$i++) {

    print "Classifying using model.$i ...\n";
    system("$SVM_CLASSIFY $svm_para_file $SVM_MODEL_DIR/model.$i $svm_pred_base.$i");
    
}
#}

my %pred=();
for(my $i=1;$i<=5;$i++) {

   
    open(FILE,"$svm_pred_base.$i");
    while(<FILE>) {
	chomp;
	push(@{$pred{$i}},$_);
    }
    close(FILE);
}




my @mean=();
my @std=();
for(my $j=0;$j<scalar(@{$pred{1}});$j++) {
    my @vec=();
    for(my $i=1;$i<=5;$i++) {
	push(@vec,$pred{$i}[$j]);

    }
    push(@mean,mean(@vec));
    push(@std,std(@vec));
    

}
my @CAs=`grep CA $pdb.orig`;
my @id=();
foreach my $ca(@CAs) {
    my $chain=substr($ca,21,1);
    my $resno=substr($ca, 22, 4);
    my $id="$resno:$chain";
    push(@id,$id);
}

my $global_quality=mean(@mean);

my $email_body="";
$email_body="Link to local quality plot: http://wallner.theophys.kth.se/output/$jobid/local.png\n";
$email_body.="Raw files: http://wallner.theophys.kth.se/output/$jobid/\n";
$email_body.="Predicted quality [0,1] , 0-worst, 1-best\n\n";
$email_body.="Global Quality: ";
$email_body.=sprintf("%-8.3f\n",$global_quality);
$email_body.="Local Quality:\n";
$email_body.=sprintf("%-5s %-8s %-8s %-8s\n","Pos","S","dist","resnum:chain");
open(OUT,">$svm_pred_base.all");
print OUT "# first column pos, second column average of five last. third column std of five last, fourth column resno:chain\n";
my $d0=3;
for(my $j=0;$j<scalar(@{$pred{1}});$j++) {
    printf OUT ("%-5d %-10.7f %-10.7f %-8s ",$j+1,$mean[$j],$std[$j],$id[$j]);
    my $m=$mean[$j];
    $m=1 if($m>1);
    $m=0 if($m<0);
    my $d=15;
    if($m>0.03846) {
	$d=sqrt(1/$m-1)*$d0;
    }

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
print OUT "set xlabel 'sequence position'\n";
print OUT "set ylabel 'predicted quality, 1-best 0-worst'\n";
print OUT "plot '$svm_pred_base.all' u 1:2 w l lw 1 title 'raw','' u 1:2 w l lw 2 smooth bezier title 'smoothed'\n";
close(OUT);
system("/opt/local/bin/gnuplot $svm_pred_base.all.gnu");
`cp $svm_pred_base.png $folder/local.png`;

if($sendemail) {
#send email with result.
    
    my $email=`cat $email_file | head`;
    my $target_name=`cat $folder/name|head`;
    
    chomp($email);
    chomp($target_name);
    if(length($target_name)>0) {
	$target_name="for $target_name";
    }
    open(EMAIL,">$folder/results.email");
    print EMAIL "From: ProQM\@cbr.su.se\n";
    print EMAIL "Subject: Local Quality Prediction ProQM-server Results $target_name\n";
    print EMAIL "To: $email\n";
    print EMAIL $email_body;
    close(EMAIL);
    
    

    `cat $folder/results.email | /usr/sbin/sendmail -t`;
} else {

    `cp  $svm_pred_base.all $pdb.ProQM`;
}

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
