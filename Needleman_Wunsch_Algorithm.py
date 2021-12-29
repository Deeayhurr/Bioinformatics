import numpy as np
import os 


class Needleman_algorithm:

    def __init__(self):
        self.gap_penalty = None
        self.same_award = None
        self.difference_penalty = None
        self.max_number_paths = None
        self.max_seq_length = None
        self.seq_counter = 0
        self.length_sequence = None
        self.length_sequence_2 = None
        
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
            
            Example:
            -2
            2
            -1
            20
            3

            The Sequence should be in Fasta format

            The input to the program should be in this format: seq1.fasta, seq2.fasta, config.txt, output.txt
        '''

        print("\t" + self.esc('31;1;4') +  "-----User Guide----"+ self.esc(0) + INFO_MESSAGE)

    def file_length(self,file_name):
        with open(file_name) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def getConfigFile(self):
        """ This function gets the input from the user and checks if everything is accurate """

        while True:
            try:
                """This is a try block for accepting the path to the config file. It also validates the if the file and contents of the file"""

                config_file = str(input(f'Input path to config file: '))

                if os.path.splitext(config_file)[1] == ".txt":
                    config_len = self.file_length(config_file)
                    if config_len == 5:
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
                                    
                            except ValueError:
                                print("The contents of the file should be integers")
                    else:
                        raise ValueError("The config file should contain 5 lines ")
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


    def getFirstSeq(self):
        """ This function gets the input from the user and checks if everything is accurate """

        while True:
            
            try:
                """This is a try block for accepting the path to the first sequence. It also validates the if the file and contents of the file"""
                
                first_sequence = str(input(f'Input path to first sequence: '))
                if  (os.path.splitext(first_sequence)[1] == ".txt"):
                    with open(first_sequence, "r") as f_seq:
                        seq_1 = f_seq.read()
                        sequences_1 = seq_1.split(">")    
                        for sequence in sequences_1:
                            position = sequence.rfind(']')
                            if position < sequence.find('\n'):
                                sequence_1 = sequence[position+1:].strip()
                                self.length_sequence = len(sequence_1.strip())
                                if self.length_sequence > self.max_seq_length:
                                    raise ValueError(f"The sequence length is {self.length_sequence} and it  has exceeded" 
                                    f"{self.max_seq_length} the maximum number of sequence length." 
                                    "Please input a sequence below or equal the maximum sequence length")
                                else:
                                    print(print(f'{self.length_sequence} is the number of charareters in {sequence_1}'))
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
                break
        return sequence_1

    def getSecondSeq(self):
        while True:  
             
            try:
                """This is a try block for accepting the path to the second sequence. It also validates the if the file and contents of the file"""

                second_sequence = str(input(f'Input path to second sequence: '))
                if (os.path.splitext(second_sequence)[1] == ".txt"):
                    with open(second_sequence, "r") as s_seq:
                        seq_2 = s_seq.read()
                        sequences_2 = seq_2.split(">")    
                        for sequence in sequences_2:
                            position = sequence.rfind(']')
                            if position < sequence.find('\n'):
                                sequence2 = sequence[position+1:].strip()
                                self.length_sequence_2 = len(sequence2.strip())
                                if self.length_sequence_2 > self.max_seq_length:
                                    raise ValueError(f"The sequence length is {self.length_sequence_2} and it  has exceeded" 
                                    f"{self.max_seq_length} the maximum number of sequence length." 
                                    "Please input a sequence below or equal the maximum sequence length")
                                else:
                                    print(print(f'{self.length_sequence_2} is the number of charareters in {sequence2}'))
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
                break
        return sequence2

    def getOutFile(self,sequence_1,sequence2):
        while True:
            try:
                """This is a try block for accepting the path to the output file. It also validates the if the file and contents of the file"""

                output_file = str(input(f'output file name: '))

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
        zero_matrix = np.zeros((self.length_sequence+1,self.length_sequence_2+1))
        return sequence_1, sequence2, zero_matrix ,output_file
        
                

    def needleman_algo(self,seq_info):
        sequence_1 = seq_info[0]
        sequence_2 = seq_info[1]
        matrix = seq_info[2]
        output_file = seq_info[3]

        res = np.empty(matrix.shape, dtype=list)

        #fill out the first row and first column of the matrix
        
        for col_index in range(1,matrix.shape[1]):
            matrix[0,col_index] =  matrix[0,col_index-1] + self.gap_penalty
            
        for row_index in range(1,matrix.shape[0]):
            matrix[row_index,0] =  matrix[row_index-1,0] + self.gap_penalty 
            
        #filling out the remaining parts of the matrix
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
        results = []

        def get_path(row, column, recorded_path = []):
            if self.seq_counter == self.max_number_paths:
                return
            if row == 0  and column == 0:
                print()
                recorded_path.append(sequence_1[0] + sequence_2[0])
                if recorded_path:
                    results.append(recorded_path)
                    self.seq_counter = self.seq_counter +1
            if res[row+1, column+1]:
                if 'D' in res[row + 1, column+1]:
                    recorded_path.append(sequence_2[column]+sequence_1[row])
                    get_path(row-1, column-1, recorded_path.copy())
                    recorded_path.pop()
                if 'L' in res[row+1, column+1]:
                    recorded_path.append(sequence_2[column]+ '-')
                    get_path(row, column-1, recorded_path.copy())
                    recorded_path.pop()
                if 'T' in res[row+1, column+1]:
                    recorded_path.append('-'+sequence_1[row])
                    get_path(row-1, column, recorded_path.copy())
                    recorded_path.pop()

        get_path(matrix.shape[0]-2,matrix.shape[1]-2)
        

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
    needleman1.getConfigFile()
    seq1 = needleman1.getFirstSeq()
    seq2 = needleman1.getSecondSeq()
    outfile = needleman1.getOutFile(seq1,seq2)
    calculated_score = needleman1.needleman_algo(outfile)