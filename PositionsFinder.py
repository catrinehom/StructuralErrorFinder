#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Program: PositionFinder
# Description: This is program is a part of the tool StructuralErrorFinder
# Version: 1.0
# Author: Catrine Høm

"""
# Import libraries
from argparse import ArgumentParser
from ErrorHandling import log, OpenFile
import random
import sys

###########################################################################
# FUNCTIONS
###########################################################################

def FindPositions(filename, min_alignment_length, logstring, outputname):
    """
    Description: Finds positions of start and stop alignment of a read which has been aligned to a contig. 
    This information is found in the blast_results.out file after using the tool BLAST.
    
    Parameters: 
        filename: Path to the blast_results.out output file from BLAST
        logstring: String to be written to log
    """
    positions = dict()
    

    file = OpenFile(filename,'r',logstring)
    pos_count = 0
    
    for line in file:
        l = line.split('\t')
        
        if len(l) == 12:
            
            alignmentlength = int(l[3])
                    
            pos1 = int(l[8])
            pos2 = int(l[9])
                   
            if alignmentlength > min_alignment_length:
                if pos1 in positions:
                    positions[pos1] += 1
                    pos_count += 1
                else:
                    positions[pos1] = 1
                    pos_count += 1
                        
                if pos2 in positions:
                    positions[pos2] += 1
                    pos_count += 1
                else:
                    positions[pos2] = 1
                    pos_count += 1
        else:
            message = 'Error! Positions is not found in this line: \n {}'.format(line)
            logstring += message
            print(message)
    
    # Get number of reads
    no_reads = pos_count//2
    message = 'The number of reads mapped to your reference is: {}\n'.format(no_reads)
    logstring += message
    print(message)   
    
    # Test if no positions were found
    if no_reads == 0:
        message = "Error! No positions found in file. \n"
        logstring += message
        log(logstring,outputname)
        sys.exit(message)
        
    return no_reads, positions, logstring


def SimulatePositions(no_reads, assembly_length):
    """
    Description: Simulates positions of start and stop alignment of a read to a contig. 
    
    Parameters:
        no_reads: Number of reads you want to simulate
        assembly_length: Length of assembly where the reads is spread out
    """
    positions = dict()
    
    for x in range (0, no_reads):
        pos1 = random.randint(0, assembly_length)
        pos2 = random.randint(0, assembly_length)
        
        if pos1 in positions:
            positions[pos1] += 1
        else:
            positions[pos1] = 1
            
        if pos2 in positions:
            positions[pos2] += 1
        else:
            positions[pos2] = 1
            
    return positions


def SumPositions(positions,stepsize):
    """
    Description: Sums position over a step<ize
    
    Parameters:
        positions: Dictionary with key as bp positions in assembly,
    and value as count of read start/stopping at that position.
    """
    summed_positions = dict()
    
    ## TODO: probably no use for this function
    ## TODO: What about first and last positions?
    
    for key in positions:
        addition = stepsize//2
        summed_positions[key+addition] = positions[key]
        for i in range(1,stepsize+1):
            if key+i in positions:
                summed_positions[key+addition] += positions[key+i]
    
    
    return summed_positions

                
def FindAssemblyLength(assemblyname, positions, logstring):
    """
    Description: Finds length of assembly
    
    Parameters:
        assemblyname: Assembly to find length of
        positions: Dictionary with count per positions in assembly
        logstring: String to be written to log
    """
    file = OpenFile(assemblyname,'r',logstring)
    
    fasta = file.read()

    try:
        header, proseq = fasta.split('\n', 1)
        seq = proseq.replace('\n', '')
        assembly_length = len(seq)
    except ValueError:
        message = 'Sequence length not found from fasta file! '\
                  'Assembly length is estimated from the read mapped at the highest position.'
        logstring += message
        print(message)
        
        assembly_length = max(positions.keys())
    
    return assembly_length, logstring


def FindDistribution(positions, logstring, printdist=False):
    """
    Description: Finds distribution of counts per position. 
    
    Parameters:
        positions: Dictionary with count per positions in assembly
        logstring: String to be written to log
        printdist: Option of print distribution to logstring, default is false.
    """
    dist = dict()
    
    for key in positions:
        count = positions[key]
        if count in dist:
            dist[count] += 1
        else:
            dist[count] = 1
    
    if printdist == True:
        message = 'This is the distribution of start/stop positions:\n'
        for key in sorted(dist):
            message += '{} {}\n'.format(key, dist[key])
        logstring += message
        print(message)  
    
    return dist, logstring
    

def FindMedian(values):
    """
    Description: Finds median of values
    
    Parameters:
        values: Values to find median from
    """
    n = len(values) 
    values.sort() 
  
    if n % 2 == 0: 
        median1 = values[n//2] 
        median2 = values[n//2 - 1] 
        median = (median1 + median2)/2
    else: 
        median = values[n//2] 
    return median


def FindOutlierValue(no_reads, assembly_length, repetitions, level, logstring):
    """
    Description: Finds outlier values of distribution based on simulation.
    
    Parameters:
        no_reads: Number of reads mapped to assembly
        assembly_length: Length of assembly
        repetitions: Number of repetitions in simulation
        logstring: String to be written to log
    """
    maxvalues = list()
    
    for i in range(0,repetitions):
        positions = SimulatePositions(no_reads, assembly_length)
        dist, logstring = FindDistribution(positions, logstring)
        maxvalues.append(max(dist.keys()))
    
    # Find outlier level
    if level == 1:
        outlier_value = max(maxvalues)
    elif level == 2:
        outlier_value = FindMedian(maxvalues)
    elif level == 3:
        outlier_value = min(maxvalues)
    
    # Print results
    message = 'Outlier value (upper limit): {}\n'.format(outlier_value)
    message += 'Positions which have higher count than the outlier value, will be seen as outliers. \n'
    logstring += message
    print(message)

    return outlier_value, logstring


def FindOutliers(positions, assembly_length, outlier_value, logstring):
    """    
    Description: Find oositions that have counts that is higher than the outlier_value.
    
    Parameters:
        positions: Dictionary with count per positions in assembly
        assembly_length: Length of assembly
        outlier_value: Upper limit for counts for them to be considered randomly distributed. 
        logstring: String to be written to log
    """
    
    # Identify outliers
    outliers = dict()
    
    # Iterate over all the items in dictionary
    for (key, value) in positions.items():
        # If item is an outlier, add to new dict
        if value > outlier_value:
            outliers[key] = value
    
    # Remove first and last position
    min_value = 1
    max_value = assembly_length
    
    if min_value in outliers:
        del outliers[min_value]
        
    if max_value in outliers:
        del outliers[max_value]    
        
    # Print results
    if bool(outliers) == True:
        message = 'This is the outliers:\n'
        
        message += 'pos\tcount\n'
        for outlier in sorted(outliers):
            message += '{}\t{}\n'.format(outlier, outliers[outlier])
    
        logstring += message
        print(message)
    
    else:
        message = 'No outlier positions found.\n'
        logstring += message
        print(message)
           
    return logstring


if __name__== "__main__":
       
###########################################################################
# GET INPUT
###########################################################################

    # Parse input from command line
    parser = ArgumentParser(description='TODO')
    parser.add_argument("-i", dest="inputname", help="input file (blast_result.out) from BLAST")
    parser.add_argument("-a", dest="assemblyname", help="assembly fasta file")
    parser.add_argument("-o", dest="outputname", help="Output filename")
    #parser.add_argument("-s", dest="stepsize", help="Step-size", type=int, default = 4)
    parser.add_argument("-m", dest="min_alignment_length", help="Minimum alignment length for reads", type=int, default = 5000)
    parser.add_argument("-l", dest="level", help="level of output (1, 2 or 3)", type=int, default = 1)
    args = parser.parse_args()

    # Define input as variables
    inputname = args.inputname
    outputname = args.outputname
    assemblyname = args.assemblyname
    #stepsize = args.stepsize
    min_alignment_length = args.min_alignment_length
    level = args.level
    
    # Define logstring
    logstring = str()
    
    # Repetitions of simulated distribution
    repetitions = 100
    
    ## Test input
    #inputname = 'blast_results.out'
    #assemblyname = 'assembly.fasta'
    #outputname = 'test_blast'

###########################################################################
# RUN PROGRAM
###########################################################################

    ### Find positions
    no_reads, positions, logstring = FindPositions(inputname, min_alignment_length, logstring, outputname) 
        
    ### Find distribution
    dist, logstring = FindDistribution(positions, logstring, printdist=True)
    
    ### Find assembly length
    assembly_length, logstring = FindAssemblyLength(assemblyname, positions, logstring)

    ### Find outlier values
    outlier_value, logstring = FindOutlierValue(no_reads, assembly_length, level, repetitions, logstring)

    ### Filter outliers
    logstring = FindOutliers(positions, assembly_length, outlier_value, logstring)
    
    ### Write to logfile
    log(logstring,outputname)
    
    
    
### PIPES IN .sh ??     
    
    ### Sum positions
    #summed_positions = SumPositions(positions, stepsize)
    #FindDistribution(summed_positions)
    #logstring = FindOutliers(summed_positions, logstring)

