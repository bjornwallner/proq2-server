#!/bin/bash
#wget http://bioinfo.ifm.liu.se/ProQ2/cache.tar.gz
wget http://bioinfo.ifm.liu.se/ProQ2/DB.tar.gz
wget http://bioinfo.ifm.liu.se/ProQ2/sspro4.db.tar.gz
tar -xvzf sspro4.db.tar.gz
rm sspro4.db.tar.gz
tar -xvzf DB.tar.gz
rm DB.tar.gz
cd apps/sspro4/
./configure.pl
cd ../..
#cd apps/svm_light/
#make
#cd ../../
