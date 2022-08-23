#!/bin/bash

for file in *.gz
do
	fd=${file%%_*}
	if [ -d $fd ]; then
		if [ -e $file ]; then
		echo "Moving file $file to folder $PWD/$fd/"
		mv $file "$PWD/$fd/"
		fi
	else
	mkdir "$PWD/$fd"
	echo "Moving file $file to folder $PWD/$fd/"
	mv $file "$PWD/$fd/"
	fi
done
