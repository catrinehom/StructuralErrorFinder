#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program: ErrorHandling
Description: This is program is a part of the tool StructuralErrorFinder
Version: 1.0
Author: Catrine Høm
"""

# Import libraries
import sys
import os
import gzip
from argparse import ArgumentParser

###########################################################################
# FUNCTIONS
###########################################################################

def log(logstring,outputname):
    """
    Description: Writes a logstring to logfile
    
    Parameters:
        logstring: String to be written to log
        outputname: Path and name to where the log is saved
    """
    # Open log file
    logname = "{}/{}.log".format(outputname, outputname)
    logfile = open(logname,"a+")
    
    # Append logstring to logfile
    logfile.write(logstring)
    logfile.close()

def CheckExistence(filenames,logstring):
    """
    Description: Check if valid filenames is provided.
    
    Parameters:
        filenames: files to be checked to exist
        logstring: string to be written to log
    """
    for file in filenames:
        if os.path.exists(file) == False:
            message = "Input Error: {} does not exist in path.".format(file)
            logstring += message
            log(logstring,outputname)
            sys.exit(message)

def CheckGZip(filename):
    """
    Description: Checks if the input file is gzipped.
    
    Parameters:
        filename: Filename to be checked if gzipped.
    """
    gzipped_type = b"\x1f\x8b"

    infile = open(filename,"rb")
    filetype = infile.read(2)
    infile.close()
    
    if filetype == gzipped_type:
        return True
    else:
        return False

def OpenFile(filename,mode,logstring):
    """
    Description: Opens the input file in chosen mode.
    
    Parameters:
        filename: file to be opened
        mode: how it should be opened
        logstring: string to be written to log
    """
    try:
        if CheckGZip(filename):
            infile = gzip.open(filename,mode)
        else:
            infile = open(filename,mode)
    except IOError as error:
        message = "Can't open file, reason: {} \n".format(error)
        logstring += message
        log(logstring,outputname)
        sys.exit(message)
    return infile

def CheckFastq(filenames,logstring):
    """
    Description: Checks if all input files (list) are in fastq format.    
    
    Parameters:
        filenames: files to be checked to contain fastq reads
        logstring: string to be written to log
    """
    fastq = list()
    fastq_type = b"@"

    # Open all files and get the first character
    for infile in filenames:
        f = OpenFile(infile, "rb",logstring)
        first_char = f.read(1);
        f.close()
        # Check if fastq
        if first_char == fastq_type:
            fastq.append(True)
        else:
            fastq.append(False)
            
    for i in range(0,len(fastq)):
        if fastq[i] == False:
            message = "Input Error: {} is a wrong format. Should be fasta format.".format(fastq[i])
            logstring += message
            log(logstring,outputname)
            sys.exit(message)

def CheckFasta(filenames,logstring):
    """
    Description: Checks if all the input files (list) are in fasta format.
    
    Parameters:
        filenames: files to be checked to contain fasta reads
        logstring: string to be written to log
    """
    fasta = list()
    fasta_type = b">"

    # Open file and get the first character
    for infile in filenames:
        f = OpenFile(infile, "rb",logstring)
        first_char = f.read(1);
        f.close()
        # Check if fasta
        if first_char == fasta_type:
            fasta.append(True)
        else:
            fasta.append(False)

    for i in range(0,len(fasta)):
        if fasta[i] == False:
            message = "Input Error: {} is a wrong format. Should be fasta format.".format(fasta[i])
            logstring += message
            log(logstring,outputname)
            sys.exit(message)

###########################################################################
# GET INPUT
###########################################################################

if __name__ == '__main__':
    # Parse input to program
    parser = ArgumentParser(description='TODO')
    parser.add_argument("-f", dest="fastqfilename", help="Fastq file", nargs = "+")
    parser.add_argument("-a", dest="assembly",help="Assembly file", nargs = "+")
    parser.add_argument("-o", dest="outputname", help="Output filename")
    parser.add_argument("-s", dest="stepsize",help="Step-size",type=int)
    parser.add_argument("-m", dest="min_alignment_length", help="Minimum alignment length for reads", type=int, default = 5000)
    parser.add_argument("-l", dest="level", help="level of output (1, 2 or 3)", type=int, default = 1)
    args = parser.parse_args()

    # Define input as variables
    fastqfilename = args.fastqfilename
    assemblyfilename = args.assembly
    outputname = args.outputname

    # Define logstring
    logstring = ''

###########################################################################
# TEST INPUT
###########################################################################

    # Test if fasta references and fastq files is exists in folder
    CheckExistence(fastqfilename,logstring)
    CheckExistence(assemblyfilename,logstring)

    # Test if assembly and input files is correct format
    CheckFasta(assemblyfilename,logstring)
    CheckFastq(fastqfilename,logstring)
            
    # Save to logfile
    log(logstring,outputname)
