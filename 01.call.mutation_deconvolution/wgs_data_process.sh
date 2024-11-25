#!/bin/bash
#

set -o nounset
set -o errexit

if [ $# -lt 1 ]	## no parameters
then
	echo "Usage: $0 [options] -o output.prefix -1 read1.fq [-2 read2.fq]" > /dev/stderr
	echo > /dev/stderr
	echo "Options:" > /dev/stderr
	echo "  -g genome  Set genome. Default hg38. Supports hg38.HBV.EBV, hg38.Hieff, hg19, mm10, rn6" > /dev/stderr
	echo "  -s size    Set minimum read size. Default: 36" > /dev/stderr
	echo "  -t thread  Set running threads. Default: auto" > /dev/stderr
	echo "  -k kit     Set kit for trimming adaptors. Default: illumina" > /dev/stderr

	exit 2
fi

samtools=/usr/bin/samtools

# default parameters
species=hg38
minSize=36
cpunum=8
seqKit=illumina
output=""
read1=""
read2=""

# read command line parameters
while getopts ":g:s:o:1:2:t:k" OPTION
do
	case $OPTION in
		g)species="$OPTARG"
			;;
		s)minSize="$OPTARG"
			;;
		t)cpunum="$OPTARG"
			;;
		o)output="$OPTARG"
			;;
		1)read1="$OPTARG"
			;;
		2)read2="$OPTARG"
			;;
		k)seqKit="$OPTARG"
			;;
		?)echo -e "\n\n***** ERROR: unsupported option detected. *****\n"
			;;
	esac
done

reference=/usr/genomes/hg38_p13/hg38p13
BWA=/usr/bin/bwa
picard=/usr/bin/picard.jar

if [ -z "$output" ]
then
	echo "Error: No output file!"
	exit 102
fi

[ -d $output ] || mkdir $output 
cd $output

echo "PWD: $PWD"   > $output.run.log
echo "CMD: $0 $@" >> $output.run.log

seqmode=""
if [ -z "$read2" ]
then
	seqmode="se"
else
	seqmode="pe"
fi

if [ $cpunum == 0 ]	## use all threads
then
	cpunum=`cat /proc/cpuinfo | grep processor | wc -l`
fi
echo -e "Species\t$species" >> $output.run.log
echo -e "minSize\t$minSize" >> $output.run.log
echo -e "seqMode\t$seqmode" >> $output.run.log
echo -e "seqKit\t$seqKit"   >> $output.run.log

echo "Running Ktrim ..."
if [ $seqmode == "pe" ]
then
	[ -s "$output.ktrim.trim.log" ] || $prgbase/ktrim -1 $read1 -2 $read2 -t $cpunum  -p $PHRED -o $output.ktrim -s $minSize -k $seqKit
	CMD="$BWA mem -t $cpunum -M -Y  -R $header $reference  $output.ktrim.read1.fq  $output.ktrim.read2.fq | $samtools view -Sb -o  $output.bam - && $samtools sort  -@ $cpunum -O bam -o  $output.sort.bam $output.bam 2> $output.bwa.log"
else
fi

echo "Running BWA in $seqmode mode ..."
echo $CMD -s $output.sort.bam ] || sh -c "$CMD"

echo "Mark duplicates ..."
	
[ -s $output.mkdup.bam ] || java -Xmx6g -Djava.io.tmpdir=./ -XX:ParallelGCThreads=$cpunum  -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $picard MarkDuplicates REMOVE_DUPLICATES=true I=$output.sort.bam O=$output.mkdup.bam METRICS_FILE=$output.mkdup.bam.mat TMP_DIR=tmp 2>$output.mkdup.log 

if [ -s $output.mkdup.bam.bai ]
then
        echo "SKIP SAM to BAM conversion"
else
	echo "samtools index ..."
	$samtools index -@ $cpunum $output.mkdup.bam 
	$samtools stats -@ $cpunum $output.mkdup.bam > $output.mkdup.bam.stats
fi
echo "Done: the main output is '$output.mkdup.bam' (with mkdup)."
