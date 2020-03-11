#!/usr/bin/env bash

# Program: StructuralErrorFinder
# Description: StructuralErrorFinder is a pipeline to find structural errors in your assembly based on your reads.
# Version: 1.0
# Author: Catrine Høm

# Usage:
    ## StructuralErrorFinder [-f <fastq reads>] [-a <fasta assembly] [-o <output name]

# This pipeline consists of 2 steps:
    ## STEP 1: BLAST READS AGAINST ASSEMBLIES
    ## STEP 2: FIND POSITIONS

###########################################################################
# GET INPUT
###########################################################################

# Start timer for logfile
SECONDS=0

# How to use program
usage() { echo "Usage: $0 [-f <fastq reads>] [-a <fasta assembly] [-o <output name]" 1>&2; exit 1; }

# Default values
stepsize=5
threads=8
min_alignment_length=5000
level=1

# Parse flags
while getopts ":f:a:o:m:l:" opt; do
    case "${opt}" in
        f)
            fastqreads=${OPTARG}
            ;;
        a)
            assembly=${OPTARG}
            ;;
        o)
            outputname=${OPTARG}
            ;;
        m)
            min_alignment_length=${OPTARG}
            ;;
        l)
            level=${OPTARG}
            ;;
        h)
            usage
            ;;
        *)
            echo "Invalid option: ${OPTARG}"
            usage
            ;;
    esac
done
shift $((OPTIND -1))

# Check if required flags are empty
if [ -z "${fastqreads}" ] || [ -z "${assembly}" ] || [ -z "${outputname}" ]; then
    echo "f, a and o are required flags."
    usage
fi

# Make output directory
[ -d $outputname ] && echo "Output directory: ${outputname} already exists. Files will be overwritten."  | tee -a $log || mkdir $outputname

# Make logfile and empty it
log=${outputname}/${outputname}.log
touch $log
cat /dev/null > $log

# Starting program
date=$(date "+%Y-%m-%d %H:%M:%S")
echo "Starting StructuralErrorFinder ($date)" | tee -a $log
echo "--------------------------------------------------" | tee -a $log
echo "StructuralErrorFinder is a pipeline to find structural errors in your assembly based on your reads." | tee -a $log
echo "" | tee -a $log

# Check format and that the files exists
python3 ErrorHandling.py -f $fastqreads -a $assembly -o $outputname -s $stepsize -m $min_alignment_length -l $level

# Check if python script exited with an error
if [ $? -eq 0 ]
then
  echo "Error handling of input done." | tee -a $log
else
  echo "Script exited due to input error." | tee -a $log
  exit 1
fi

# Print files used
echo "Input fastq files are ${fastqreads}" | tee -a $log
echo "Assembly used is ${assembly}" | tee -a $log

echo "Time stamp: $SECONDS seconds." | tee -a $log
echo "" | tee -a $log

### Convert fastq reads into fasta format
if [[ $fastqreads =~ \.gz$ ]];
then
    gunzip $fastqreads
    fastqreads="${fastqreads%.*}"
    sed -n '1~4s/^@/>/p;2~4p' $fastqreads > reads.fasta
    fastareads=reads.fasta
    gzip $fastqreads
else
    sed -n '1~4s/^@/>/p;2~4p' $fastqreads > reads.fasta
    fastareads=reads.fasta
fi

###########################################################################
# STEP 1:  BLAST READS AGAINST ASSEMBLIES
###########################################################################

echo "Starting STEP 1: BLAST reads against assembly." | tee -a $log

[ -d $outputname/databases ] && echo "Output directory: ${outputname} already exists. Files will be overwritten."  | tee -a $log || mkdir $outputname/database

# Make database of assembly and BLAST reads against the database
db=./databases/assembly_db
blast_result="blast_results.out"

makeblastdb -in $assembly -dbtype nucl -out $db
blastn -db $db -query $fastareads -out $blast_result -outfmt 6

echo "Time stamp: $SECONDS seconds." | tee -a $log
echo "" | tee -a $log

###########################################################################
# STEP 2: FIND POSITIONS
###########################################################################

echo "Starting STEP 3: Find positions" | tee -a $log
python3 PositionFinder.py -i $blast_result -a $assembly -o $outputname -m $min_alignment_length -l $level

# Check if python script exited with an error
if [ $? -eq 0 ]
then
  echo "Positions with high probability of structural errors found." | tee -a $log
else
  echo "Script exited due to error." | tee -a $log
  exit 1
fi

echo "Time stamp: $SECONDS seconds." | tee -a $log
