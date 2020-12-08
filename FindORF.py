############################################
# Create a command line program that will  #
# find all ORFs in a fasta file with 1     #
# or many sequences.                       #
# input: A fasta file with DNA sequence(s).#
# output: A tab file with predicted ORFs.  #
############################################

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from FindORF_utils import find_orfs

descript = """
            This is a program that will find all ORF(s) in a FASTA file with 1 or many sequences. It does so by using the 33 tables from NCBI. The following url will give you             more information about each table: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

           """

parser=argparse.ArgumentParser(description=descript)

parser.add_argument( 
    "-i",
    '--inputfile', 
    type=str, 
    metavar="<fasta file>", 
    help="Enter the DNA sequence(s) in FASTA format",
    required = True
)
parser.add_argument(
    "-t",
    '--table', 
    type=int, 
    metavar="<Integer>", 
    help="Table which will be used to translate DNA sequences",
    required = True
)
parser.add_argument(
    "-l",
    '--length', 
    type=int, 
    metavar="<Integer>", 
    help="The minimum length of the ORF Protein Sequences",
    required = True
)


    
args = vars(parser.parse_args())

with open("results.tsv", "w") as file:
    for fasta_sequence in SeqIO.parse(args['inputfile'], 'fasta'):
        name = fasta_sequence.id
        sequence = fasta_sequence.seq
        orfs = find_orfs(sequence, args['table'], args['length'])
        file.write("Seq\tStrand\tFrame\tP_len\n")
        for i in orfs:
            file.write(f"{i.protein}\t{i.strand}\t{i.frame}\t{i.length}\n")

