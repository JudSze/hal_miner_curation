#!/bin/bash

enzyme=$1
input_dir=$2
test_file=$3

cd $input_dir
echo $PWD

for FILE in $PWD/$enzyme/crossval_fasta/*.fasta
do
    mafft $FILE > ${FILE}_aligned.afa
    hmmbuild ${FILE}.hmm ${FILE}_aligned.afa 
done

mv $PWD/${enzyme}/crossval_fasta/*.hmm $PWD/${enzyme}/crossval_profile

for FILE in $PWD/${enzyme}/crossval_profile/*.hmm
do
    echo "hmmsearch running"
    hmmsearch $FILE $test_file > $FILE.out
    echo "$FILE.out"
done

mv $PWD/${enzyme}/crossval_profile/*.out $PWD/${enzyme}/crossval_hmmsearch_res