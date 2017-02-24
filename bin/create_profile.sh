#!/bin/bash
#module add blast/2.0.11
# This is a simple script which will carry out all of the basic steps
# required to make a PSIPRED V2 prediction. Note that it assumes that the
# following programs are in the appropriate directories:
# blastpgp - PSIBLAST executable (from NCBI toolkit)
# makemat - IMPALA utility (from NCBI toolkit)
# psipred - PSIPRED V2 program
# psipass2 - PSIPRED V2 program

# The name of the BLAST data bank
#set dbname = /afs/pdc.kth.se/home/a/arnee/structpredict/SCOP/scop-1.57/data/nrdb90+scop
full_script_path=`readlink -f $0`
script_dir=`dirname $full_script_path`
dir=`dirname $script_dir`
dbname=$dir/DB/uniref90.fasta


# Where the NCBI programs have been installed


#set ncbidir = /Users/bjornw/Research/Apps/blast-2.2.16/bin/
ncbidir=$dir/apps/blast-2.2.18_x86_64/bin/
if [ ! -f ".ncbirc" ] 
    then
    echo "[ncbi]" > .ncbirc
    echo "Data=$dir/apps/blast-2.2.18_x86_64/data/" >> .ncbirc
fi
echo $dbname

#set ncbidir = /usr/ebiotools/bin/



# Where the PSIPRED V2 programs have been installed
execdir=$dir/apps/psipred25/bin/

# Where the PSIPRED V2 data files have been installed
datadir=$dir/apps/psipred25/data/

basename="${1%.*}"
outfile=$basename.mtx

if  [ -f  $outfile ] 
then
    echo Already exist $outfile ...
#exit
#else
#echo "Let's go!!!!"
fi

echo "OUT: $outfile"
echo "BASE: $basename"
#exit;


#\cp -f $1 psitmp.fasta
echo cp -f $1 psitmp.$$.fasta
cp -f $1 psitmp.$$.fasta

echo "Running PSI-BLAST with sequence" $1 "..."
echo "$ncbidir/blastpgp -a 4 -j 3 -h 0.001 -d $dbname -F F -i psitmp.$$.fasta -C psitmp.$$.chk -Q psitmp.$$.psi"

$ncbidir/blastpgp -a 4 -j 3 -h 0.001 -d $dbname -F F -i psitmp.$$.fasta -C psitmp.$$.chk -o psitmp.$$.blastpgp -Q psitmp.$$.psi > /dev/null #& $rootname.blast


echo "Running Makemat..."
echo psitmp.$$.chk > psitmp.$$.pn
echo psitmp.$$.fasta > psitmp.$$.sn
$ncbidir/makemat -P psitmp.$$
#cp psitmp.$$.mtx $rootname.mtx
#cp psitmp.$$.chk $rootname.chk
#cp psitmp.$$.psi $rootname.psi
#echo Cleaning up ...
#rm -f psitmp.$$.* error.log 
#exit;
echo "Predicting secondary structure..."
echo Pass1 ...

$execdir/psipred psitmp.$$.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > psitmp.$$.ss

echo Pass2q ...

$execdir/psipass2 $datadir/weights_p2.dat 1 1.0 1.0 psitmp.$$.ss2 psitmp.$$.ss 
cp psitmp.$$.ss2 $basename.ss2
cp psitmp.$$.mtx $basename.mtx
cp psitmp.$$.chk $basename.chk
cp psitmp.$$.psi $basename.psi
cp psitmp.$$.blastpgp $basename.fasta.blastpgp # this is for accpro will skip the psiblast runs
# Remove temporary files

echo Cleaning up ...
rm -f psitmp.$$.* error.log
exit;

#rm $rootname.ss $rootname.blast $rootname.ss2 #$rootname.horiz

#echo "Final output files:" $rootname.ss2 $rootname.horiz
echo "Final output files:" $rootname.horiz $rootname.ss2 
echo "Finished."
