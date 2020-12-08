############################################
# Create a command line program that will  #
# find all ORFs in a fasta file with 1     #
# or many sequences.                       #
# input: A fasta file with DNA sequence(s).#
# output: A fasta file with predicted ORFs.#
############################################
from Bio import SeqIO

class ORF:
    def __init__(self, protein, strand, frame):
        self.protein = protein
        self.strand = strand
        self.frame = frame
        self.length = len(protein)
        
#find ORFs
def find_orfs(sequence, minimum_length, table):
    orfs = []
    for strand, nucleotides in [(+1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):
            length = 3 * (len(sequence)-frame) 
            for proteins in nucleotides[frame:frame+length].translate(table).split("*"):
                if len(proteins) >= minimum_length:
                    if(strand == +1):
                        orfs.append(ORF(proteins, "+", frame+1))
                    elif(strand ==-1):
                        orfs.append(ORF(proteins, "+", frame+1))
                            
                                 
    return orfs
