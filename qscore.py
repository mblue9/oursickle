# Input:    Bio.SeqIO iterator
# Output:   A sequence object with class fastq_record extended from Bio.SeqIO
#   Create class Fastq_Record
#   Add phred_quality
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class Fastq_Record(SeqRecord):
    '''A class to extend a SeqRecord to contain the human readable QScores
    and add function to trim the reads
    '''
    def __init__(self, record):
        self.__dict__.update(record.__dict__)
        self.Qscore = self.letter_annotations["phred_quality"]

    def trim(self, five_prime, three_prime):
        '''A function to trim a read
        takes as input two integers >= 0 and extracts and trims the sequence to these indices
        '''        
        self = self[five_prime:three_prime]


