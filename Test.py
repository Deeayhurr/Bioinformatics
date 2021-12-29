import pytest
from Needleman_Wunsch_Algorithm import Needleman_algorithm
import numpy as np


def test_needleman_algo():
    """
    TESTING needelman algo using mock data
    """

    needleman_algorithm = Needleman_algorithm()
    needleman_algorithm.gap_penalty = -2
    needleman_algorithm.same_award = 2
    needleman_algorithm.difference_penalty = -1
    needleman_algorithm.max_number_paths = 3

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