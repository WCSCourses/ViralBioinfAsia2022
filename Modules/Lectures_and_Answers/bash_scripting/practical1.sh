#!/bin/bash

while read line || [ -n "$line" ]; #This will make sure that the file's last line is read by the while loop.
do
	line_array=($line)
	name=${line_array[0]}
	seq=${line_array[1]}
	echo -e ">$name \n${seq:8}"
done < $1
