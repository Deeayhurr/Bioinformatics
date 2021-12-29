import numpy as np
import os 
from Bio import SeqIO
from Bio import Entrez


class Needleman_algorithm:

    def __init__(self):
        self.gap_penalty = None
        self.same_award = None
        self.difference_penalty = None
        self.max_number_paths = None
        self.max_seq_length = None
        self.type_of_algorithm = None
        self.seq_counter = 0
    
    def fetch_from_bio(self, protein):
        Entrez.email = "a.adeyemi13@gmail.com"
        #handle = Entrez.efetch(db="nucleotide", id="EU490707", rettype="gb", retmode="text")
        handle = Entrez.efetch(db="protein", id=protein, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        print(record.description)
        print(record.id)
        print(len(record.features))
        rec_sequence = record.seq
        print(rec_sequence)
        print(f"The length of sequence is {len(rec_sequence)}")
        return rec_sequence

    def esc(self,code):
        return f'\033[{code}m'

    def userguide(self):
        """ This is a function that has enclosed the userguide of this program """

        INFO_MESSAGE = '''

            * Needleman-Wunsch Algorithm 
            The config file shoud contain these details in the order described below:
            - GAP˙PENALTY (int)
            - SAME˙AWARD (int)
            - DIFFERENCE˙PENALTY (int)
        
            - MAX˙SEQ˙LENGTH (int)
            - MAX˙NUMBER˙PATHS (int)
            - ALGORITHM˙TYPE (int) : 1 for Needleman-Wunsch algorithm
                                     2 for Waterman-Smith algorithm
            
            Example:
            -2
            2
            -1
            20
            3
            2

            The Sequence should be in Fasta format

            The input to the program should be in this format: seq1.fasta, seq2.fasta, config.txt, output.txt
        '''

        print("\t" + self.esc('31;1;4') +  "-----User Guide----"+ self.esc(0) + INFO_MESSAGE)

    def file_length(self,file_name):
        with open(file_name) as f:
            for i, l in enumerate(f):
                pass
        return i + 1
    def config_file(self):
        while True:
            try:
                """This is a try block for accepting the path to the config file. It also validates the file and contents of the file"""

                #config_file = str(input(f'Input path to config file: '))
                config_file = "config.txt"

                if os.path.splitext(config_file)[1] == ".txt":
                    config_len = self.file_length(config_file)
                    if config_len == 6:
                        with open(config_file, "r") as conf:
                            rows = conf.readlines()
                            try:
                                for count,row in enumerate(rows) :
                                    if count== 0:
                                        self.gap_penalty = int(row)
                                        print(f"\ngap penalty is {self.gap_penalty}")
                                    if count == 1:
                                        self.same_award = int(row)
                                        print(f"same award is {self.same_award}")
                                    if count == 2:
                                        self.difference_penalty = int(row)
                                        print(f"difference penalty is {self.difference_penalty}")
                                    if count == 3:
                                        self.max_seq_length = int(row)
                                        print(f"Maximum number of sequence length is {self.max_seq_length}")
                                    if count == 4:
                                        self.max_number_paths = int(row)
                                        print(f"max_number_paths is {self.max_number_paths}")
                                    if count == 5:
                                        self.type_of_algorithm = int(row)
                                        if self.type_of_algorithm == 1:
                                            algorithm = " Needleman-Wunsch algorithm"
                                        elif self.type_of_algorithm == 2:
                                            algorithm = "Waterman-Smith algorithm"
                                            print(f"Type of algorithm is {algorithm}")
                                        else:
                                            raise ValueError("The 6th line should be either 1 or 2 for  Needleman-Wunsch and Waterman-Smith algorithm respectively  ")          
                            except ValueError:
                                print("The contents of the file should be integers")
                    else:
                        raise ValueError("The config file should contain 6 lines ")
                else:
                    raise ValueError("The config file should be in the .txt file format ")
            except ValueError as err:
                print(err.args)
                continue
            except UnboundLocalError:
                print("The file is empty")
                continue
            except FileNotFoundError:
                print("Thie file does not exist")
                continue
            else :
                break
    
    def get_first_sequence(self):
        """This is a method for accepting the path to the first sequence. It also validates the if the file and contents of the file"""
        while True:

            try:
                #first_sequence = str(input(f'Input path to first sequence: '))
                #first_sequence = "seq_1.txt"
                first_sequence = "sequence.txt"
                if  (os.path.splitext(first_sequence)[1] == ".txt"):
                    with open(first_sequence, "r") as f_seq:
                        seq_1 = f_seq.read()
                        sequences_1 = seq_1.split(">")    
                        for sequence in sequences_1:
                            position = sequence.rfind(']')
                            if position < sequence.find('\n'):
                                sequence_1 = sequence[position+1:].strip()
                                length_sequence = len(sequence_1.strip())
                                if length_sequence > self.max_seq_length:
                                    raise ValueError(f"The sequence length is {length_sequence} and it  has exceeded" 
                                    f"{self.max_seq_length} the maximum number of sequence length." 
                                    "Please input a sequence below or equal the maximum sequence length")
                                else:
                                    print(print(f'{length_sequence} is the number of charareters in {sequence_1}'))
                else:
                    raise ValueError("The extension for the first sequence should be in .txt ")
            except ValueError as err:
                print(err.args)
                continue
            except UnboundLocalError:
                print("The file is empty")
                continue
            except FileNotFoundError:
                print("Thie file does not exist")
                continue
            else :
                return sequence_1

    def get_second_sequence(self):
        while True:  
             
            try:
                """This is a try block for accepting the path to the second sequence. It also validates the if the file and contents of the file"""

                #second_sequence = str(input(f'Input path to second sequence: '))
                second_sequence =  "sequence-3.txt"
                if (os.path.splitext(second_sequence)[1] == ".txt"):
                    with open(second_sequence, "r") as s_seq:
                        seq_2 = s_seq.read()
                        sequences_2 = seq_2.split(">")    
                        for sequence in sequences_2:
                            position = sequence.rfind(']')
                            if position < sequence.find('\n'):
                                sequence2 = sequence[position+1:].strip()
                                length_sequence_2 = len(sequence2.strip())
                                if length_sequence_2 > self.max_seq_length:
                                    raise ValueError(f"The sequence length is {length_sequence_2} and it  has exceeded" 
                                    f"{self.max_seq_length} the maximum number of sequence length." 
                                    "Please input a sequence below or equal the maximum sequence length")
                                else:
                                    print(print(f'{length_sequence_2} is the number of charareters in {sequence2}'))
                else:
                    raise ValueError("The extension for the second sequence should be in .txt")
            except ValueError as err:
                print(err.args)
                continue
            except UnboundLocalError:
                print("The file is empty")
                continue
            except FileNotFoundError:
                print("Thie file does not exist")
                continue
            else :
                return sequence2

    def getSequenceinfo(self):
        """ This function gets the input from the user and checks if everything is accurate """
        
        self.config_file() 
       
        while True:
            try:
                read_from_bio = str(input(f'Do you wanna read from text or from biopython reply "yes" to read from .txt and "no" to read from biopython: '))
                read_from_bio.lower
                if read_from_bio == "no":
                    #protein_1  =  str(input(f'input protein ID:'))
                    protein_1 ="Q4FZY2"
                    sequence_1 = self.fetch_from_bio(protein_1) 
                    #protein_2  =  str(input(f'input protein ID:'))
                    protein_2 = "Q8N6T7"
                    sequence2 = self.fetch_from_bio(protein_2) 
                elif read_from_bio == "yes":
                    sequence_1 = self.get_first_sequence()
                    sequence2 = self.get_second_sequence()
                else:
                    raise ValueError(f'The expected input is "yes" of "no" ')
            except ValueError as err:
                print(err.args)
                continue
            else:
                break
        
        while True:
            try:
                """This is a try block for accepting the path to the output file. It also validates the if the file and contents of the file"""

                #output_file = str(input(f'output file name: '))
                output_file = "o_wsp.txt"

                if os.path.splitext(output_file)[1] == ".txt":
                    break
                else:
                    raise ValueError("The output file name should be in the .txt file format ")
            except ValueError as err:
                print(err.args)
                continue
            else :
                break
        #Generate the matrix and fill it out with zeroes
        zero_matrix = np.zeros((len(sequence_1)+1,len(sequence2)+1))
        return sequence_1, sequence2, zero_matrix ,output_file
                
    
    def needleman_algo(self,seq_info):
        sequence_1 = np.asarray(seq_info[0])
        sequence_2 = np.asarray(seq_info[1])


        matrix = seq_info[2]
        output_file = seq_info[3]

        res = np.empty(matrix.shape, dtype=list)

        #fill out the first row and first column of the matrix
        if self.type_of_algorithm == 1: 
            for col_index in range(1,matrix.shape[1]):
                matrix[0,col_index] =  matrix[0,col_index-1] + self.gap_penalty
                
            for row_index in range(1,matrix.shape[0]):
                matrix[row_index,0] =  matrix[row_index-1,0] + self.gap_penalty
        else:
            for col_index in range(1,matrix.shape[1]):
                matrix[0,col_index] = 0
                
            for row_index in range(1,matrix.shape[0]):
                matrix[row_index,0] = 0

            
        #filling out the remaining parts of the matrix for needleman-wunsch algorithm
        if self.type_of_algorithm == 1:
            for row_index in range(1,matrix.shape[0]):
                for col_index in range(1,matrix.shape[1]):
                    if sequence_1[row_index-1] == sequence_2[col_index-1]:
                        diagonal_score = matrix[row_index - 1][col_index - 1] + self.same_award
                    else:
                        diagonal_score = matrix[row_index - 1][col_index - 1] + self.difference_penalty
                    top_score = matrix[row_index - 1][col_index] + self.gap_penalty
                    left_score = matrix[row_index][col_index - 1] + self.gap_penalty

                    max_val = max(diagonal_score, top_score, left_score)
                    matrix[row_index][col_index] = max_val
            
                    temp_res = []
                    if diagonal_score == max_val:
                        temp_res.append('D')
                    if top_score == max_val:
                        temp_res.append('T')
                    if left_score == max_val:
                        temp_res.append('L')

                    res[row_index][col_index] = temp_res
            score = matrix[-1][-1] 
            print(matrix)
        else:
            #filling out the remaining parts of the matrix for waterman-smith algorithm
            for row_index in range(1,matrix.shape[0]):
                for col_index in range(1,matrix.shape[1]):
                    if sequence_1[row_index-1] == sequence_2[col_index-1]:
                        diagonal_score = matrix[row_index - 1][col_index - 1] + self.same_award
                    else:
                        diagonal_score = matrix[row_index - 1][col_index - 1] + self.difference_penalty
                    top_score = matrix[row_index - 1][col_index] + self.gap_penalty
                    left_score = matrix[row_index][col_index - 1] + self.gap_penalty

                    if diagonal_score < 0:
                        diagonal_score = 0 
                    if top_score < 0:
                        top_score = 0 
                    if left_score < 0:
                        left_score = 0

                    max_val = max(diagonal_score, top_score, left_score)
                    matrix[row_index][col_index] = max_val
            
                    
                    temp_res = []
                    if diagonal_score == max_val:
                        temp_res.append('D')
                    if top_score == max_val:
                        temp_res.append('T')
                    if left_score == max_val:
                        temp_res.append('L')

                    res[row_index][col_index] = temp_res
            score = np.amax(matrix)
            print(matrix)

        #The function for getting the needleman-wunsch algorithm path
        results = []
        def get_needleman_path(row, column, recorded_path = []):
            if self.seq_counter == self.max_number_paths:
                return
            if row == 0  and column == 0:
                recorded_path.append(sequence_1[0] + sequence_2[0])
                if recorded_path: 
                    results.append(recorded_path)
                    self.seq_counter = self.seq_counter +1
                        
            if res[row+1, column+1]:
                if 'D' in res[row + 1, column+1]:
                    recorded_path.append(sequence_1[row]+sequence_2[column])
                    get_needleman_path(row-1, column-1, recorded_path.copy())
                    recorded_path.pop()
                if 'L' in res[row+1, column+1]:
                    recorded_path.append( '-'+sequence_2[column])
                    get_needleman_path(row, column-1, recorded_path.copy())
                    recorded_path.pop()
                if 'T' in res[row+1, column+1]:
                    recorded_path.append(sequence_1[row]+'-')
                    get_needleman_path(row-1, column, recorded_path.copy())
                    recorded_path.pop()

        #The method for getting the waterman-smith algorithm path               
        def get_Waterman_path(row, column, recorded_path = []):
            maxElement = np.amax(matrix)
            if self.seq_counter == self.max_number_paths:
                return
            if matrix[row][column] == 0:
                recorded_path.pop()
                #recorded_path.append(sequence_1[row-1] + sequence_2[column-1])
                if recorded_path: 
                    results.append(recorded_path)
                    self.seq_counter = self.seq_counter +1
            if matrix[row][column] == maxElement:

                recorded_path.append(sequence_1[row-1]+sequence_2[column-1])    
                if res[row, column]:
                    if 'D' in res[row , column]:
                        recorded_path.append(sequence_1[row-2]+sequence_2[column-2])
                        get_Waterman_path(row-1, column-1, recorded_path.copy())
                        recorded_path.pop()
                    if 'L' in res[row, column]:
                        recorded_path.append('-'+sequence_2[column-2])
                        get_Waterman_path(row, column-1, recorded_path.copy())
                        recorded_path.pop()
                    if 'T' in res[row, column]:
                        recorded_path.append(sequence_1[row-2]+'-')
                        get_Waterman_path(row-1, column, recorded_path.copy())
                        recorded_path.pop()
            if res[row, column]:
                if 'D' in res[row , column]:
                    recorded_path.append(sequence_1[row-2]+sequence_2[column-2])
                    get_Waterman_path(row-1, column-1, recorded_path.copy())
                    recorded_path.pop()
                if 'L' in res[row, column]:
                    recorded_path.append('-'+sequence_2[column-2])
                    get_Waterman_path(row, column-1, recorded_path.copy())
                    recorded_path.pop()
                if 'T' in res[row, column]:
                    recorded_path.append(sequence_1[row-2]+'-')
                    get_Waterman_path(row-1, column, recorded_path.copy())
                    recorded_path.pop()
    
        if self.type_of_algorithm == 1:
            #Getting the path for the needleman-wunsch algorithm
            get_needleman_path(matrix.shape[0]-2,matrix.shape[1]-2)
        else:
            #Getting the path for the waterman-smith algorithm
            max_element_index = np.where(matrix == np.amax(matrix))
            max_element_index = list(zip(max_element_index[0], max_element_index[1]))
            for index_list in range(len(max_element_index)):
                row = max_element_index[index_list][0]
                column = max_element_index[index_list][1]
                get_Waterman_path(row,column)


            
            

        with open (f'{output_file}', 'w+') as f:
            f.write(f'Score : {score} \n\n')
            for count, result in enumerate(results):
                if count < self.max_number_paths:
                    f.write('\n-----\n')
                    seq_1 = [val[0] for val in reversed(result)]
                    seq_2 = [val[1] for val in reversed(result)]
                    f.writelines(''.join(seq_1)+'\n')
                    f.writelines(''.join(seq_2)+'\n')
        

# The main program
if __name__ == '__main__':
    needleman1 = Needleman_algorithm()
    needleman1.userguide()
    Algo_info = needleman1.getSequenceinfo()
    calculated_score = needleman1.needleman_algo(Algo_info)
    

    


