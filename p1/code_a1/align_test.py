"""
Unit tests provide a way of testing individual components of your program.

This will help you with debugging and making sure you don't break your code!


Here, we provide some unit tests and a skeleton for a few others.
Note that you do not have to follow these exactly, they are designed to help you.


What other unit tests might you want to write?
 - Think about the traceback and writing the output file. 
 - Try to write at least one or two additional tests as you think of them.


To run:
  python align_test.py

Make sure align.py is located in the same directory, and the test_example.input file is present!
"""

import unittest

from align import *

TEST_INPUT_FILE="test_example.input"

class TestAlignmentClasses(unittest.TestCase):

    def test_match_matrix(self):
        """
        Tests match matrix object
        """
        match_matrix = MatchMatrix()
        match_matrix.set_score("A", "C", 5)
        self.assertEqual(match_matrix.get_score("A", "C"), 5)

    def test_score_matrix_score(self):
        """
        Tests score matrix object score set + get methods
        """
        score_matrix = ScoreMatrix('M', 3, 3)
        score_matrix.set_score(2, 1, 5)
        self.assertEqual(score_matrix.get_score(2, 1), 5)

    def test_score_matrix_pointers(self):
        """
        Tests score matrix object pointer set + get methods
        """
        score_matrix = ScoreMatrix('M', 2, 2)
        score_matrix.set_pointers(1, 1, 0, 0, 'M')
        score_matrix.set_pointers(1, 1, 0, 0, 'Ix')
        self.assertEqual(score_matrix.get_pointers(1, 1), [(0,0,'M'), (0,0,'Ix')])

    def test_param_loading(self):
        """
        Tests AlignmentParameters "load_params_from_file()" function
        """
        align_params = AlignmentParameters()
        align_params.load_params_from_file(TEST_INPUT_FILE)
        self.assertEqual(align_params.seq_a, "AATGC")
        self.assertEqual(align_params.seq_b, "AGGC")
        self.assertTrue(align_params.global_alignment)
        self.assertEqual(align_params.dx, 0.1)
        self.assertEqual(align_params.ex, 0.5)
        self.assertEqual(align_params.dy, 0.6)
        self.assertEqual(align_params.ey, 0.3)
        self.assertEqual(align_params.alphabet_a, "ATGC")
        self.assertEqual(align_params.alphabet_b, "ATGCX")
        self.assertEqual(align_params.len_alphabet_a, 4)
        self.assertEqual(align_params.len_alphabet_b, 5)

        # test that the match match is set up correctly
        #  if this fails, make sure you are loading the asymmetric matrix properly!
        match_mat = align_params.match_matrix
        self.assertEqual(match_mat.get_score("A", "X"), 0.3)
        self.assertEqual(match_mat.get_score("C", "G"), -0.3)
        self.assertEqual(match_mat.get_score("G", "C"), 0)


    def test_update_ix(self):
        """
        Test AlignmentAlgorithm's update Ix
        """
        

        # configure alignment params
        align_params = AlignmentParameters()
        align_params.dy = 1
        align_params.ey = 0.5

        # create an alignment object
        align = Align("", "")
        align.align_params = align_params
        align.m_matrix = ScoreMatrix("M", 5, 4)
        align.ix_matrix = ScoreMatrix("Ix", 5, 4)
        align.m_matrix.set_score(2,2, 3)
        align.ix_matrix.set_score(2,2, 2.5)

        # run the method!
        pointers = align.ix_matrix.get_pointers(3, 2)
        align.update_ix(3, 2)

        score = align.ix_matrix.get_score(3,2)
        self.assertEqual(score, 2)
        ### FILL IN for checking pointers! ###
        # note - in this example, it should point to M -AND- Ix
        # check by hand!
        self.assertEqual(align.ix_matrix.get_pointers(3, 2), [(2,2,'M'), (2,2,'Ix')])


    def test_update_m(self):
        """
        Test AlignmentAlgorithm's update M
        """
        ### FILL IN ###
        return
    
    def test_update_iy(self):
        """
        Test AlignmentAlgorithm's update Iy
        """
        # configure alignment params
        align_params = AlignmentParameters()
        align_params.dx = 1
        align_params.ex = 0.5

        # create an alignment object
        align = Align("", "")
        align.align_params = align_params
        align.m_matrix = ScoreMatrix("M", 5, 4)
        align.iy_matrix = ScoreMatrix("Iy", 5, 4)
        align.m_matrix.set_score(2,2, 3)
        align.iy_matrix.set_score(2,2, 2.5)

        # run the method!
        pointers = align.iy_matrix.get_pointers(2, 3)
        align.update_iy(2, 3)
        

        score = align.iy_matrix.get_score(2,3)
        self.assertEqual(score, 2)
        ### FILL IN for checking pointers! ###
        # note - in this example, it should point to M -AND- Ix
        # check by hand!
        self.assertEqual(align.iy_matrix.get_pointers(2, 3), [(2,2,'M'), (2,2,'Iy')])

    def test_traceback_start(self):
        """
        Tests that the traceback finds the correct start
        Should test local and global alignment!
        """
        # configure alignment params
        align_params = AlignmentParameters()
        align_params.dy = 1
        align_params.ey = 0.5
        align_params.len_alphabet_a = 4
        align_params.len_alphabet_b = 4
        align_params.seq_a = "ATGC"
        align_params.seq_b = "ATCG"

        align_params.global_alignment = True

        # create an alignment object
        align = Align("", "")
        align.align_params = align_params
        align.m_matrix = ScoreMatrix('M', 5, 5)
        align.ix_matrix = ScoreMatrix('Ix', 5, 5)
        align.iy_matrix = ScoreMatrix('Iy', 5, 5)

        for i in np.arange(1, 5):
            for j in np.arange(1, 5):
                align.m_matrix.set_score(i, j, i*j)
                align.ix_matrix.set_score(i, j, 0)
                align.iy_matrix.set_score(i, j, i*j)

        self.assertEqual(align.find_traceback_start(), (16, [(4, 4)]))


        for i in np.arange(1, 5):
            for j in np.arange(1, 5):
                align.m_matrix.set_score(i, j, i*j)
                align.ix_matrix.set_score(i, j, 50)
                align.iy_matrix.set_score(i, j, i*j)
        self.assertEqual(align.find_traceback_start(), (50, [(4, 1), (4, 2), (4, 3), (4, 4), (1, 4), (2, 4), (3, 4)]))

        align_params.global_alignment = False
        for i in np.arange(1, 5):
            for j in np.arange(1, 5):
                align.m_matrix.set_score(i, j, 1/(1 + i*j))
        self.assertEqual(align.find_traceback_start(), (0.5, [(1,1)]))


if __name__=='__main__':
    unittest.main()