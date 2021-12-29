import numpy as np
import os 
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq as sq


class Multiple_Sequence_Allignmnet:

    def fetch_file(self):
       with open("msa_Sequences.txt", "w+") as handle:
        for record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
            sequence = record.seq
            reverse_complement = sequence.reverse_complement()
            transcribe_seq = reverse_complement.transcribe()
            translated_seq =  transcribe_seq.translate()
            handle.write('\n-----\n')
            handle.writelines(''.join((translated_seq).replace("*",""))+'\n')

        

            












# The main program
if __name__ == '__main__':
    msa = Multiple_Sequence_Allignmnet()
    Algo_info = msa.fetch_file()
    