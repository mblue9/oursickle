#!/usr/bin/python

import sys, getopt
from Bio import SeqIO
from qscore import Fastq_Record
from window_size import WindowSize
from window import slide_window

# exception for when invalid Fastq format specified
class InvalidFastqformatException(Exception):
    pass

# main function
def main(argv):
    # get variables
    fastq_format = ''
   
    # use getopt to get options
    try:
        opts, args = getopt.getopt(argv,"hi:f:m:o:t:",["input_fastq", "fastq_format=", "min_len=", "output_fastq", "threshold="])
    except:
        print 'sickle.py -f <[solexa|illumina|sanger]> -m <min_len> -t <threshold> input_fastq output_fastq'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
           print 'sickle.py -f <[solexa|illumina|sanger]> -m <min_len> -t <threshold> input_fastq output_fastq'
           sys.exit()
        elif opt in ("-f", "--fastq_format"):
           fastq_format = arg
        elif opt in ("-m", "--min_len"):
           min_len = int(arg)
        elif opt in ("-t", "--threshold"):
           threshold = int(arg)
        elif opt in ("-i", "--input_fastq"):
           input_fastq = arg
        elif opt in ("-o", "--output_fastq"):
           output_fastq = arg
    
    
    # mapping from fastq_format to Bio.Seq
    types = {"sanger": 'fastq', "illumina": 'fastq-illumina', "solexa": 'fastq-solexa'}
    
    # check min_len is valid integer > 0
    if not isinstance(min_len, int) or min_len <= 0:
        raise Exception("The specified min_len must be an integer > 0")
    
    # check fastq_format is valid
    if fastq_format not in types:
        raise InvalidFastqformatException("You have specified an invalid Fastq format: " + fastq_format + ". Available formats: solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (CASAVA >= 1.8).")

    # check threshold is valid integer >= 1 and <= 40
    if not isinstance(min_len, int) or min_len < 1 or min_len > 40:
        raise Exception("The specified threshold must be an integer >= 1 and <= 40")

    # open the file, Bio.SeqIO will handle exceptions
    handle = open(input_fastq, "rU")
    
    # open file to write to, Bio.SeqIO will handle exceptions
    out_handle = open(output_fastq, "w")
    
    for record in SeqIO.parse(handle, types[fastq_format]):
        # determine window size
        window_size = WindowSize(len(record.seq))
        
        # determine cut points
        five_prime, three_prime = slide_window(record.letter_annotations["phred_quality"], window_size, threshold)
        
        # if whole read poor quality indexes will both be -1 
        if five_prime == -1 and three_prime == -1:
            continue
        
        # trim the reads
        record = record[five_prime:three_prime]
        
        # output read
        SeqIO.write(record, out_handle, "fastq")
    
if __name__ == "__main__":
   main(sys.argv[1:])
