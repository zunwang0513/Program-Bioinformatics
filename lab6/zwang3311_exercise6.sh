#!/bin/bash


#To set these variables to 0 means that in default no realignment of the bam file, the output vcf will not be gunzipped, the verbose mode is turned off, and there will be no index of the output bam file. The default output name is study.vcf

realign=0
gunzip=0
v=0
index=0
answer=0
output="study.vcf"
helpmessage="This is a pipeline to have the function of Mapping, Alignment Improvement, and Variant Calling. There are some options you can use to modify this pipeline. You will need to input two input reads, a reference file, a millsFile, and the output filename. The output will be a vcf file containing the variants and two txt: SNP.txt and indels.txt that contains the SNP and indel information.\n\t
	-a File location for read1\n\t
	-b File location for read2\n\t
	-r File location for reference genome\n\t
	-e Perform Read Realignment\n\t
	-o the name of the output file\n\t
	-f File location for mills file\n\t
	-z output VCF file will be gunzipped\n\t
	-v Turn on verbose mode\n\t
	-i Index the output BAM file\n\t
	-h To print this help message
	"
#################################GetOpts#######################################
while getopts a:b:r:eo:f:zvih option
do
	case $option in
		a) reads1=$OPTARG;;
		b) reads2=$OPTARG;;
		r) ref=$OPTARG;;
		e) realign=1;;
		o) output=$OPTARG;;
		f) millsFile=$OPTARG;;
		z) gunzip=1;;
		v) v=1;;
		i) index=1;;
		h) echo -e $helpmessage; exit;;
	esac
done

##############################FILE CHECK#######################################

#Check if read1 and read2 exist. This is to check the existance of read1 and read2 at the same time.

if [ ! -f $reads1 ]
then
	echo "ERROR: reads1 does not exist"
	if [ ! -f $reads2 ]
       	then
		echo "ERROR: reads2 does not exist"
	fi
	echo "Program Exit"
	exit
fi

if [ ! -f $reads2 ]
then
	echo "ERROR: read2 does not exist"
	echo "Program Exit"
	exit
fi

#Check if reference file exist.

if [ ! -f $ref ]
then
	echo "ERROR: reference file does not exist"
	echo "Program Exit"
	exit
fi

#Check if there is already a file named as the output filename. I write this while loop to ensure that if the user input something else than y or n, the warning message will pop up until the user input y or n.

if [ -f $output ]
then

while [ $answer != 'y' ]
do
	echo "Warning: Output VCF File already exists. Permission to overwrite [y/n]?"
	read answer
	if [ $answer == 'n' ]; then
		exit
	fi
done
fi

#Check if millsFIle exists.

if [ ! -f $millsFile ]
then
	echo "ERROR: Mills File does not exist."
	exit
fi


#These files must exist to proceed, a,b,r,and f flags must have input
if [ -z $reads1 ]
then
        echo "ERROR: no input of read1"
	exit
fi
if [ -z "$reads2" ]
then
        echo "ERROR: no input of read2"
	exit
fi
if [ -z "$ref" ]
then
        echo "ERROR: no input of ref"
	exit
fi
if [ -z "$millsFile" ]
then
        echo "ERROR: no input of millsFile"
	exit
fi
#################################Mapping#######################################

#Index the reference file using bwa.

if [ v ]
then
	echo "Indexing reference"
fi

bwa index $ref

if [ v ]
then
	echo "Index reference complete"
fi

#Mapping the reads to a file called lane.sam and add header.

if [ v ]
then
        echo "Aligning the reads to the reference"
fi

bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' $ref $reads1 $reads2 > lane.sam

if [ v ]
then
        echo "Alignment complete"
fi

#Clean up read pairing info and flags

if [ v ]
then
        echo "Cleaning up the read pairing info and flags"
fi


samtools fixmate -O bam lane.sam lane_fixmate.bam

if [ v ]
then
        echo "Cleaning up complete"
fi


#Sort the bam file from name order into coordinate order

if [ v ]
then
        echo "Sorting bam file from name order into coordinate order"
fi


samtools sort -O bam -o lane_sorted.bam -T /tmp/lane_temp lane_fixmate.bam

if [ v ]
then
        echo "Sorting complete"
fi


#To prepare for alignment improvement, we need to create fai and dict file for the reference. Also, the bam file needed to be indexed. So if the user need Improvement of the alignment. User should use index function.

if [ v ]
then
        echo "Creating fai file for reference"
fi


samtools faidx $ref

if [ v ]
then
        echo "Creating dict file for reference"
fi

#java -jar picard-2.10.0_picard.jar CreateSequenceDictionary R=$ref
samtools dict $ref -o ${ref%.*}.dict

if [ $index ]
then

	if [ v ]
	then
        	echo "indexing the bam file"
	fi

	samtools index lane_sorted.bam
fi

################################Improvement####################################

#Realignment process of the sorted bam file.

if [ $realign ]
then
	if [ v ]
	then
        	echo "Realigning the bam file"
	fi

	java -Xmx2g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I lane_sorted.bam -o lane.intervals --known $millsFile 2> zunwang0513.log 

	java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I lane_sorted.bam -targetIntervals lane.intervals -known $millsFile -o lane_realigned.bam 2>> zunwang0513.log
	if [ v ]
	then
        	echo "Realignment complete"
	fi

fi

###############Other Improvement Not Neccesary########################
#java -Xmx4g -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R $ref -knownSites knownSite -I lane_realigned.bam -o lane_recal.table
#java -Xmx2g -jar GenomeAnalysisTK.jar -T PrintReads -R $ref -I lane_realigned.bam --BSQR lane_recal.table -o lane_recal.bam
###########################################################################

#index the file for variant calling

if [ $index ]
then

	if [ v ]
	then
        	echo "Indexing the bam file"
	fi

	samtools index lane_realigned.bam
fi

#############################Calling Variant###################################

#Use bcftools to mpileup and call the variant. There are realigned or gunzip option that makes slightly difference.

if [ v ]
then
        echo "Calling variants from the bam file"
fi


if [ $realign ]
then
	bcftools mpileup -Ob -o study.bcf  -f $ref lane_realigned.bam

else
	bcftools mpileup -Ob -o study.bcf  -f $ref lane_sorted.bam
fi

bcftools call -vmO z -o $output.gz study.bcf

if [ $gunzip == 0 ]
then
	gunzip $output.gz
fi

if [ v ]
then
        echo "Variant Calling Complete. Your file is ready"
fi



###############################NOT NEEDED DELETED#############################

##########################Generate SNP and indel txt###########################

# The following is to generate bed file from vcf file and split the SNP and indels into two txt files.

#if [ $gunzip ]
#then
#	gunzip -c $output.gz > $output
#fi

#sed '/^#/d' $output | cut -f1,2,4,5 > temp.vcf

#awk -vOFS='\t' '{$5=$2+length($4)-length($3);$6=length($4)-length($3)}1' temp.vcf | cut -f1,2,5,6 > temp2.vcf

#awk '{gsub("chr", "");print}' temp2.vcf > temp3.vcf

#awk '$4 == 0 {print > ("SNP1.txt"); next} {print > ("indels1.txt")}'  temp3.vcf

#rm temp.vcf temp2.vcf temp3.vcf

#if [ $gunzip ]
#then
#	rm $output
#fi






