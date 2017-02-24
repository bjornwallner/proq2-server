#!/bin/sh
if [ $# -ne 2 ]
then
	echo "need two parameters:seq_file, output_file." 
	exit 1
fi
/tmp/leffe/proq2-server/apps/sspro4/script/generate_flatblast.pl /tmp/leffe/proq2-server/apps/sspro4/blast2.2.8/ /tmp/leffe/proq2-server/apps/sspro4/script/ /tmp/leffe/proq2-server/apps/sspro4/data/big/big_98_X /tmp/leffe/proq2-server/apps/sspro4/data/nr/nr $1 $2 
