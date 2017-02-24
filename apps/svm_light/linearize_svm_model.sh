#!/bin/bash

model=$1
modelout=$model.linear


/afs/pdc.kth.se/home/b/bjornw/CASP9/apps/svm_light/svm_classify /dev/null $model /dev/null | grep LINEAR_WEIGHTS > /scratch/svm.weights.$$
head -n 11 $model |sed s/[0-9]\ #\ kernel\ type/4\ #\ kernel\ type/ > /scratch/svm.head.$$

cat /scratch/svm.head.$$ /scratch/svm.weights.$$ > $modelout

echo "Output written to $modelout"
rm -f /scratch/svm.head.$$
rm -f /scratch/svm.weights.$$