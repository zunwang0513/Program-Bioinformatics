#!/bin/bash

filenum=10
seqnum=8
verbose=0
#By default 10 files will be created and 8 sequences will be created, variable verbose is to keep track of whether v is used.

###############################GETOPTS#########################################
while getopts n:m:v option
do
	case $option in 
		n)filenum=$OPTARG;;
		m)seqnum=$OPTARG;;
		v)verbose=1;;
	esac
done

#########################Delete Previous FIles#################################
#Previously I know there are exactly 10 files, so I only need to delete exactly 10 files. Now the number of existing and new created files are dynamic, so I need to delete all the fasta files begin with 'seq'.
if [ $verbose ]
then
	find . -type f -name 'seq*.fasta' -print -delete | sort -n | xargs -l echo 'deleting'

#the verbose option, to print the deleting file

else
	find . -type f -name 'seq*.fasta' -delete
fi

###########################Create New Files####################################

for FILE in $(seq 1 $filenum)
do
	#if v is used, this if condition will be called to print the file processing
	if [ $verbose ]
	then
		echo "Creating seq"$FILE".fasta"
	fi

	for seq in $(seq 1 $seqnum)
	do
		echo '>seq'$FILE'_'$seq >>seq$FILE.fasta
		cat /dev/urandom | tr -dc 'ACGT' | fold -w 50 | head -n1 >> seq$FILE.fasta
		#The verbose option, to print the seq number
		if [ $verbose ]
		then
			echo "Generating seq"$seq
		fi
	done
done

