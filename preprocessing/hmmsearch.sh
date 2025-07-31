#!/bin/bash


profile=$1
input_dir=$2

cd $input_dir
echo $PWD

for FILE in $PWD/*.fasta
do
    echo "hmmsearch running"
    hmmsearch $profile $FILE > $FILE.out
    echo "$FILE.out"
done