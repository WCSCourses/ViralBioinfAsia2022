#!/bin/bash

fd=$1
ref_file=$2

files=($(ls $fd/*.gz))
file1=${files[0]}
file2=${files[1]}

## 1. Run Trimmometic ##
echo "1. Run Trimmometic: $fd"
f1_trim=${file1/.fastq.gz/.trim.fastq.gz}
f2_trim=${file2/.fastq.gz/.trim.fastq.gz}
f1_unpair=${file1/.fastq.gz/.unpair.fastq.gz}
f2_unpair=${file2/.fastq.gz/.unpair.fastq.gz}
trim_cmd="trimmomatic PE -phred33 $file1 $file2 $f1_trim $f1_unpair $f2_trim $f2_unpair LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:40"

echo "Input files: $file1 $file2"
echo -e "Output files: \n$f1_trim \n$f2_trim \n$f1_unpair \n$f2_unpair"
##----Run trimmometic command
$trim_cmd

echo -e "\n\n"

## 2. Run Minimap2
echo "2. Run Minimap2: $fd"
out_sam=${file1/_*.fastq.gz/.sam}
map_cmd="minimap2 -ax sr -o $out_sam $ref_file $f1_trim $f2_trim"

echo -e "Input files: \n$f1_trim \n$f2_trim"
echo -e "Output files: \n$out_sam"
##----Run minimap2 command
$map_cmd

echo -e "\n\n"

## 3. Run samtools
echo "3. Run samtools: $fd"
echo "-----Convert SAM to BAM---------"
out_bam=${file1/_*.fastq.gz/.bam}
echo "Input files: $out_sam"
echo "Output files: $out_bam"
bam_cmd="samtools view -Shb -o $out_bam $out_sam"
##---Run command: Convert SAM to BAM
$bam_cmd

echo -e "\n" 
echo "-----Sort BAM file-------"
sorted_bam=${file1/_*.fastq.gz/.sorted.bam}
echo "Input files: $out_bam"
echo "Output files: $sorted_bam"
sort_cmd="samtools sort -@ 2 -o $sorted_bam $out_bam"
##---Run command: Sort BAM file
$sort_cmd

echo -e "\n"
echo "-----Filter paired mapped-------"
pair_bam=${file1/_*.fastq.gz/.sorted.pair.bam}
echo "Input files: $sorted_bam"
echo "Output files: $pair_bam"
pair_cmd="samtools view -hb -f 2 -o $pair_bam $sorted_bam"
##-----Run filter paired mapped
$pair_cmd

echo -e "\n\n"

## 4. Run Flagstat
echo "4. Run Flagstat: $fd"
echo "Input files: $pair_bam"
flag_cmd="samtools flagstat $pair_bam"
##-----Run FLAGSTAT
$flag_cmd
