import pytest
from Needleman_Wunsch_Algorithm import Needleman_algorithm
import numpy as np
from io import StringIO

conf_file = StringIO('config.txt\n')

def test_config_file(monkeypatch):
    """
    TESTING needelman algo using mock data
    """

    needleman_algorithm = Needleman_algorithm()

    monkeypatch.setattr('sys.stdin', conf_file)
    needleman_algorithm.getConfigFile()

    assert needleman_algorithm.gap_penalty == -2
    assert needleman_algorithm.same_award == 2
    assert needleman_algorithm.difference_penalty == -1
    assert needleman_algorithm.max_seq_length == 550
    assert needleman_algorithm.max_number_paths == 5

first_sequence = StringIO('seq_1.txt\n')

def test_first_seq_file(monkeypatch):
    """
    TESTING needelman algo using mock data
    """
    needleman_algorithm = Needleman_algorithm()

    needleman_algorithm.gap_penalty = -2
    needleman_algorithm.same_award = 2
    needleman_algorithm.difference_penalty = -1
    needleman_algorithm.max_number_paths = 5
    needleman_algorithm.max_seq_length  = 550

    monkeypatch.setattr('sys.stdin', first_sequence)
    out_put = needleman_algorithm.getFirstSeq()

    expected = 'ATTACA'
    assert out_put == expected

sec_sequence = StringIO('seq_2.txt\n')

def test_second_seq_file(monkeypatch):
    """
    TESTING needelman algo using mock data
    """

    needleman_algorithm = Needleman_algorithm()

    needleman_algorithm.gap_penalty = -2
    needleman_algorithm.same_award = 2
    needleman_algorithm.difference_penalty = -1
    needleman_algorithm.max_number_paths = 5
    needleman_algorithm.max_seq_length  = 550

    monkeypatch.setattr('sys.stdin', sec_sequence)
    out_put = needleman_algorithm.getSecondSeq()

    expected  = "ATGCT"

   

    assert out_put == expected


out_file = StringIO('output.txt\n')
def test_out_file(monkeypatch):
    """
    TESTING needelman algo using mock data
    """
    seq1 ="ATTACA",
    seq2 = "ATGCT",
    needleman_algorithm = Needleman_algorithm()
    needleman_algorithm.length_sequence  = 6
    needleman_algorithm.length_sequence_2  = 5



    monkeypatch.setattr('sys.stdin', out_file)
    out_put = needleman_algorithm.getOutFile(seq1,seq2) 
 
    
    expected  = (
    ("ATTACA",),
    ("ATGCT",),
    np.zeros((7,6)),
    "output.txt")

    assert out_put[0] == expected[0]
    assert out_put[1] == expected[1]
    assert out_put[2].all() == expected[2].all()
    assert out_put[3] == expected[3]

def test_needleman_algo():
    """
    TESTING needelman algo using mock data
    """

    needleman_algorithm = Needleman_algorithm()
    needleman_algorithm.gap_penalty = -2
    needleman_algorithm.same_award = 2
    needleman_algorithm.difference_penalty = -1
    needleman_algorithm.max_number_paths = 5
    needleman_algorithm.max_seq_length  = 550

    seq_info = [
    "ATTACA",
    "ATGCT",
    np.zeros((7,6)),
    "output.txt"]

    with open("out.txt") as f1:
        expected = f1.read()
    
    needleman_algorithm.needleman_algo(seq_info)

    with open("output.txt") as f2:
        out_put = f2.read()

    assert out_put == expected

    