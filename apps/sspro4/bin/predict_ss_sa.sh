#!/bin/sh
#predict ss, sa for a single sequence from scratch.
if [ $# -ne 2 ]
then
	echo "need 2 parameters:seq_file(in fasta format), output_file." 
	exit 1
fi
#output: a file with predicted ss and sa.
/local/www/services/ProQ2/apps/sspro4/script/predict_ss_sa.pl /local/www/services/ProQ2/apps/sspro4/blast2.2.8/ /local/www/services/ProQ2/apps/sspro4/data/big/big_98_X /local/www/services/ProQ2/apps/sspro4/data/nr/nr /local/www/services/ProQ2/apps/sspro4/server/predict_seq_ss.sh /local/www/services/ProQ2/apps/sspro4/server/predict_seq_sa.sh /local/www/services/ProQ2/apps/sspro4/script/ $1 $2 
