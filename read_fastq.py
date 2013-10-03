#module to read in fastq file
#input is fastq file 

#from webpage
from Bio import SeqIO

def read_fastq():
record_iterator = SeqIO.parse("Quality/example.fastq", "fastq")
out_handle = open("Quality/temp.qual", "w")
SeqIO.write(record_iterator, out_handle, "qual")
out_handle.close()

#test
firstline = file.readline()
if firstline[0] != ["@"]
    raise Exception("file is not in fastq format expect @ as first character!")
