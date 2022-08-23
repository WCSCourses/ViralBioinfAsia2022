#!/bin/bash

ref_file=$1
gene_file=$2

ref=$(grep -v ">" $ref_file | tr -d "\n")

while read line || [ -n "$line" ]; #This will make sure that the file's last line is read by the while loop.
do
	line=$(echo $line | tr -d "\r")
	line_array=($line)
	gene=${line_array[0]}
	start=${line_array[1]}
	end=${line_array[2]}
	position=$((start-1))
	length=$((end-start+1))
	gene_seq=${ref:position:length}
	echo -e ">$gene \n$gene_seq"
done < $gene_file
