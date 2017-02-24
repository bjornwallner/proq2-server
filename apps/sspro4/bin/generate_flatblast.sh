#!/bin/sh
if [ $# -ne 2 ]
then
	echo "need two parameters:seq_file, output_file." 
	exit 1
fi
/local/www/services/ProQ2/apps/sspro4/script/generate_flatblast.pl /local/www/services/ProQ2/apps/sspro4/blast2.2.8/ /local/www/services/ProQ2/apps/sspro4/script/ /local/www/services/ProQ2/apps/sspro4/data/big/big_98_X /local/www/services/ProQ2/apps/sspro4/data/nr/nr $1 $2 
