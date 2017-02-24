#!/bin/bash

for f in 1 2 3 4 5; do 
    ./svm_classify svm.in /scratch/DB/ProQ2/data/model.$f leffe > /scratch/weights.$f
    head -n 11 /scratch/DB/ProQ2/data/model.$f> /scratch/head.$f
    cat /scratch/head.$f /scratch/weights.$f  > ../ProQ2/data/model.linear.$f
    rm -fr /scratch/weights.$f
    rm -fr /scratch/head.$f
done
