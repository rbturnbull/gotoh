# -*- coding: utf-8 -*-
import unittest
from gotoh import msa
import numpy as np


class GotohMSATest(unittest.TestCase):
    def setUp(self):
        pass

    # def test_null_seqi( self ):
    #     array =  np.asarray( [1,2,3] )
    #     array = np.expand_dims(array, axis=1)

    #     alignment = msa(array,None, gap_open=0.0)

    def test_seq_to_seq(self):
        x = np.asarray([1, 2, 3])
        x = np.expand_dims(x, axis=1)

        y = np.asarray([1, 2, 4])
        y = np.expand_dims(y, axis=1)

        gold = [
            [1, 1],
            [2, 2],
            [-1, 4],
            [3, -1],
        ]

        alignment = msa(x, y, gap_open=0.0)

        np.testing.assert_array_equal(alignment, gold)

    def test_seq_to_seq_flip(self):
        x = np.asarray([1, 2, 3])
        x = np.expand_dims(x, axis=1)

        y = np.asarray([1, 2, 4, 3])
        y = np.expand_dims(y, axis=1)

        gold = [
            [1, 1],
            [2, 2],
            [-1, 4],
            [3, 3],
        ]

        alignment = msa(x, y, gap_open=0.0)

        np.testing.assert_array_equal(alignment, gold)

    def test_alignment_to_seq(self):
        x = np.asarray([[1, 1], [2, 5], [4, 3]])
        y = np.asarray([1, 2, 4])
        y = np.expand_dims(y, axis=1)

        gold = [
            [1, 1, 1],
            [2, 5, 2],
            [4, 3, 4],
        ]

        alignment = msa(x, y, gap_open=0.0)
        np.testing.assert_array_equal(alignment, gold)

    def test_alignment_to_alignment(self):
        x = np.asarray([[1, 1], [2, 5], [4, 3]])
        y = np.asarray([[1, 8], [2, 5], [8, 8]])

        gold = [
            [ 1, 1,  1,  8],
            [ 2, 5,  2,  5],
            [-1,-1,  8,  8],
            [ 4, 3, -1, -1],
        ]

        alignment = msa(x, y, gap_open=-1e-8)

        np.testing.assert_array_equal(alignment, gold)


if __name__ == "__main__":
    unittest.main()
