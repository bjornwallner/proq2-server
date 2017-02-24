#!/bin/sh
#predict ss, sa for a single sequence from scratch.
if [ $# -ne 2 ]
then
	echo "need 2 parameters:seq_file(in fasta format), output_file." 
	exit 1
fi
#output: a file with predicted ss and sa.
/tmp/leffe/proq2-server/apps/sspro4/script/predict_ssa.pl /tmp/leffe/proq2-server/apps/sspro4/blast2.2.8/ /tmp/leffe/proq2-server/apps/sspro4/data/big/big_98_X /tmp/leffe/proq2-server/apps/sspro4/data/nr/nr /tmp/leffe/proq2-server/apps/sspro4/server/predict_seq_ss.sh /tmp/leffe/proq2-server/apps/sspro4/script/ $1 $2 
