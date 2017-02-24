#!/usr/bin/perl -w
use File::Basename;
$args=join(' ',@ARGV);
$infile=$ARGV[0];
$dir=dirname($infile);
open(SBATCH,">$infile.sbatch");
print SBATCH "#!/bin/bash\n";
#print SBATCH "/local/www/services/ProQ2/apps/ProQ2/bin/ProQ2.pl $args > $infile.log 2>&1 \n";
print SBATCH "/local/www/services/ProQ2/apps/ProQ2/bin/ProQ2.pl $args\n"; #> $infile.log 2>&1 \n";
close(SBATCH);
system("chmod a+x $infile.sbatch");
#$cmd="cd $dir;sbatch -v -N1 --exclusive $infile.sbatch";
#system("echo $cmd > $dir/log");
system("which sbatch >> $dir/log");
system("cd $dir;/usr/bin/sbatch -v -N1 --exclusive $infile.sbatch >> $infile.sbatch.out");

exit;
while(1)
{
    my $ProQ2_running=`ps aux|grep 'ProQ2.pl' |grep -vc grep`;

    print "ProQ2_running $ProQ2_running\n";
    chomp($ProQ2_running);
    if(not($ProQ2_running)) {
	last;
    } else {
	sleep(10);
    }
}
print "/local/www/services/ProQ2/apps/ProQ2/bin/ProQ2.pl $args\n";
system("/local/www/services/ProQ2/apps/ProQ2/bin/ProQ2.pl $args > $infile.log 2>&1 ");
