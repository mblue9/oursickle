# Input:    Bio.SeqIO iterator
# Output:   A sequence object with class fastq_record extended from Bio.SeqIO
#   Create class Fastq_Record
#   Add phred_quality
from Bio import SeqIO

class Fastq_Record(SeqRecord):
    '''A class to extend a SeqRecord to contain the human readable QScores
    '''
    Qscore = self.letter_annotations["phred_quality"]

