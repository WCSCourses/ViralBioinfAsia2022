#https://gist.github.com/astatham/621901
#https://www.biostars.org/p/105388/#

while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
	#if the line begins with > use the seq name as the filename
        outfile=${line#>}.fa
	#these seqs have / in their name so replace with underscore
	outfile=${outfile//\//_}
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < ${1}
