
"""

This file provides skeleton code for align.py. 

Locations with "FILL IN" in comments are where you need to add code.

Note - you do not need to follow this set up! It is just a suggestion, and may help for program design and testing.


Usage: python align.py input_file output_file

"""


import sys
import numpy as np


#### ------ USEFUL FUNCTIONS ------- ####
def fuzzy_equals(a, b):
    """
    Checks if two floating point numbers are equivalent.
    """
    epsilon = 10**(-6) 
    return (abs(a - b) < epsilon)
    

#### ------- CLASSES ------- ####

class MatchMatrix(object):
    """
    Match matrix class stores the scores of matches in a data structure
    """

    def __init__(self, nrow = 1, ncol = 1, alphabet_a = '', alphabet_b = ''):
        self.nrow = nrow
        self.ncol = ncol
        self.alphabet_a = alphabet_a
        self.alphabet_b = alphabet_b
        self.match_matrix = np.empty((nrow, ncol))

    def set_score(self, a, b, score):
        """
        Updates or adds a score for a specified match

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
           score = the score to set it for
        """
        if a not in self.alphabet_a:
            if len(self.alphabet_a) != 0:
                self.nrow += 1
                self.match_matrix = np.vstack(self.match_matrix, np.zeros(self.ncol))
            self.alphabet_a += a
        if b not in self.alphabet_b:
            if len(self.alphabet_b) != 0:
                self.ncol += 1
                self.match_matrix = np.hstack(self.match_matrix, np.zeros(self.nrow))
            self.alphabet_b += b
        index_a = self.alphabet_a.index(a)
        index_b = self.alphabet_b.index(b)
        self.match_matrix[index_a, index_b] = score


    def get_score(self, a, b):
        """
        Returns the score for a particular match, where a is the
        character from sequence a and b is from sequence b.

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
        Returns:
           the score of that match
        """
        index_a = self.alphabet_a.index(a)
        index_b = self.alphabet_b.index(b)
        return self.match_matrix[index_a, index_b]



class ScoreMatrix(object):
    """
    Object to store a score matrix, which generated during the alignment process. The score matrix consists of a 2-D array of
    ScoreEntries that are updated during alignment and used to output the maximum alignment.
    """
    class ScoreEntries(object):
        def __init__(self, score, pointer_list):
            self.score = score
            self.pointer_list = pointer_list

    def __init__(self, name = '', nrow = 1, ncol = 1):
        self.name = name # identifier for the score matrix - Ix, Iy, or M
        self.nrow = nrow
        self.ncol = ncol

        arr = np.array([self.ScoreEntries(None, []) if i!=0 and j!=0 \
            else self.ScoreEntries(0, []) for i in range(nrow) for j in range(ncol)])
        self.score_matrix = arr.reshape((nrow, ncol))


    def get_score(self, row, col):
        """
        Returns the score for a cell in a matrix
        Inputs:
           row = row in matrix
           col = col in matrix
        Returns:
            the score in cell (row, col) in the matrix
        """
        return self.score_matrix[row, col].score
        
    def set_score(self, row, col, score):
        """
        Sets the score for a cell in a matrix
        Inputs:
           row = row in matrix
           col = col in matrix
           score = the score in cell (row, col) in the matrix
        Updates:
            the score for cell (row, col) in the matrix
        """
        self.score_matrix[row, col].score = score

    def get_pointers(self, row, col):
        """
        Returns the indices of the entries that are pointed to
        This should be formatted as a list of tuples:
         ex. [(1,1), (1,0)]
        """
        return self.score_matrix[row, col].pointer_list

    def set_pointers(self, row, col, pointer_row, pointer_col, matrix):
        """
        Sets the pointers for a cell in a matrix
        Inputs:
           row = row in matrix
           col = col in matrix
           pointer_row = row the pointer in cell (row, col) is pointing to 
           pointer_col = col the pointer in cell (row, col) is pointing to
           matrix = matrix the pointer is pointing to
        Updates:
            the pointer list for cell (row, col) in the matrix
        """
        self.score_matrix[row, col].pointer_list.append((pointer_row, pointer_col, matrix))

    def print_scores(self):
        """
        Returns a nicely formatted string containing the scores in the score matrix. Use this for debugging!

        Example:
        M=
            0.0, 0.0, 0.0, 0.0, 0.0
            0.0, 1.0, 0.0, 0.0, 0.0
            0.0, 1.0, 1.0, 1.0, 1.0
            0.0, 0.0, 1.0, 1.0, 1.0
            0.0, 0.0, 2.0, 2.0, 1.0
            0.0, 0.0, 1.0, 2.0, 3.0

        """
        print(self.name)
        for i in np.arange(self.nrow):
            for j in np.arange(self.ncol):
                print(self.get_score(i, j), sep = ', ', end = " ", flush = True)
            print('\n')
            
    def print_pointers(self):
        """
        Returns a nicely formatted string containing the pointers for each entry in the score matrix. Use this for debugging!
        """
        print("Pointer List for matrix {}: ".format(self.name))
        print('\n')
        for i in np.arange(self.nrow):
            for j in np.arange(self.ncol):
                print("Pointer list for position [{}][{}] : {}". format(i, j, \
                self.get_pointers(i, j)))
                print('\n')
        

class AlignmentParameters(object):
    """
    Object to hold a set of alignment parameters from an input file.
    """
    def __init__(self):
        # default values for variables that are filled in by reading
        # the input alignment file
        self.seq_a = ""
        self.seq_b = ""
        self.global_alignment = True 
        self.dx = 0
        self.ex = 0
        self.dy = 0
        self.ey = 0
        self.alphabet_a = "" 
        self.alphabet_b = ""
        self.len_alphabet_a = 0
        self.len_alphabet_b = 0
        self.match_matrix = MatchMatrix(0, 0, '', '')

    def load_params_from_file(self, input_file): 
        """
        Reads the parameters from an input file and stores in the object

        Input:
           input_file = specially formatted alignment input file
        """
        with open(input_file, 'r') as params:
            self.seq_a = params.readline().strip()
            self.seq_b = params.readline().strip()
            if int(params.readline().strip()) == 1:
                self.global_alignment = False
            self.dx, self.ex, self.dy, self.ey = [float(i) for i in params.readline().strip().split()]
            self.len_alphabet_a = int(params.readline().strip())
            self.alphabet_a = params.readline().strip()
            self.len_alphabet_b = int(params.readline().strip())
            self.alphabet_b = params.readline().strip()
            self.match_matrix = MatchMatrix(self.len_alphabet_a, self.len_alphabet_b, self.alphabet_a, self.alphabet_b)
            while True:
                line = params.readline().strip().split()
                if line == []:
                    break;
                self.match_matrix.set_score(line[2], line[3], float(line[4]))
            if not self.global_alignment:
                assert np.any(self.match_matrix.match_matrix < 0)
        params.close()


class Align(object):
    """
    Object to hold and run an alignment; running is accomplished by using "align()"
    """

    def __init__(self, input_file, output_file):
        """
        Input:
            input_file = file with the input for running an alignment
            output_file = file to write the output alignments to
        """
        self.input_file = input_file
        self.output_file = output_file
        self.align_params = AlignmentParameters()
        self.align_score = 0
        self.alignments = []

        self.m_matrix = ScoreMatrix('M', self.align_params.len_alphabet_a + 1,
            self.align_params.len_alphabet_b + 1)
        self.ix_matrix = ScoreMatrix('Ix', self.align_params.len_alphabet_a + 1,
            self.align_params.len_alphabet_b + 1)
        self.iy_matrix = ScoreMatrix('Iy', self.align_params.len_alphabet_a + 1,
            self.align_params.len_alphabet_b + 1)

    def align(self):
        """
        Main method for running alignment.
        """
        # load the alignment parameters into the align_params object
        self.align_params.load_params_from_file(self.input_file)

        # populate the score matrices based on the input parameters
        self.populate_score_matrices()

        # perform a traceback and write the output to an output file
        self.align_score, self.alignments = self.traceback()
        self.write_output()

    def populate_score_matrices(self):
        """
        Method to populate the score matrices based on the data in align_params.
        Should call update(i,j) for each entry in the score matrices
        """
        self.m_matrix = ScoreMatrix('M', len(self.align_params.seq_a) + 1,
            len(self.align_params.seq_b) + 1)
        self.ix_matrix = ScoreMatrix('Ix', len(self.align_params.seq_a) + 1,
            len(self.align_params.seq_b) + 1)
        self.iy_matrix = ScoreMatrix('Iy', len(self.align_params.seq_a) + 1,
            len(self.align_params.seq_b) + 1)

        nrow = len(self.align_params.seq_a) + 1
        ncol = len(self.align_params.seq_b) + 1
        for i in np.arange(1, nrow):
            for j in np.arange(1, ncol):
                self.update(i, j)

    def update(self, row, col):
        """
        Method to update the matrices at a given row and column index.

        Input:
           row = the row index to update
           col = the column index to update
        """
        self.update_m(row, col)
        self.update_ix(row, col)
        self.update_iy(row, col)

    def update_m(self, row, col):
        """
        Method to update the M matrix at a given row and column index given 
        the alignment update rules

        Input:
           row = the row index to update
           col = the column index to update

        Updates:
            the score for cell (row, col) in the M matrix
        """
        letter_a = self.align_params.seq_a[row - 1]
        letter_b = self.align_params.seq_b[col - 1]

        m = self.m_matrix.get_score(row - 1, col - 1) + \
            self.align_params.match_matrix.get_score(letter_a, letter_b)
        ix = self.ix_matrix.get_score(row - 1, col - 1) + \
            self.align_params.match_matrix.get_score(letter_a, letter_b)
        iy = self.iy_matrix.get_score(row - 1, col - 1) + \
            self.align_params.match_matrix.get_score(letter_a, letter_b)

        max_score = max(m, ix, iy)
        if not self.align_params.global_alignment:
            max_score = max(0, max_score)
        self.m_matrix.set_score(row, col, max_score)
        
        pointer_row = row - 1
        pointer_col = col - 1
        if fuzzy_equals(max_score, m):
            self.m_matrix.set_pointers(row, col, pointer_row, pointer_col, 'M')
        if fuzzy_equals(max_score, ix):
            self.m_matrix.set_pointers(row, col, pointer_row, pointer_col, 'Ix')
        if fuzzy_equals(max_score, iy):
            self.m_matrix.set_pointers(row, col, pointer_row, pointer_col, 'Iy') 

    def update_ix(self, row, col):
        """
        Method to update the Ix matrix at a given row and column index given 
        the alignment update rules

        Input:
           row = the row index to update
           col = the column index to update

        Updates:
            the score for cell (row, col) in the Ix matrix
        """
        m = self.m_matrix.get_score(row - 1, col) - \
            self.align_params.dy
        ix = self.ix_matrix.get_score(row - 1, col) - \
            self.align_params.ey
        
        max_score = max(m, ix)
        if self.align_params.global_alignment == False:
            max_score = max(0, max_score)
        self.ix_matrix.set_score(row, col, max_score)
        
        pointer_row = row - 1
        pointer_col = col
        if fuzzy_equals(max_score, m):
            self.ix_matrix.set_pointers(row, col, pointer_row, pointer_col, 'M')
        if fuzzy_equals(max_score, ix):
            self.ix_matrix.set_pointers(row, col, pointer_row, pointer_col, 'Ix')

    def update_iy(self, row, col):
        """
        Method to update the Iy matrix at a given row and column index given 
        the alignment update rules

        Input:
           row = the row index to update
           col = the column index to update

        Updates:
            the score for cell (row, col) in the Iy matrix
        """
        m = self.m_matrix.get_score(row, col - 1) - \
            self.align_params.dx
        iy = self.iy_matrix.get_score(row, col - 1) - \
            self.align_params.ex
        
        max_score = max(m, iy)
        if self.align_params.global_alignment == False:
            max_score = max(0, max_score)
        self.iy_matrix.set_score(row, col, max_score)
        
        pointer_row = row 
        pointer_col = col - 1
        if fuzzy_equals(max_score, m):
            self.iy_matrix.set_pointers(row, col, pointer_row, pointer_col, 'M')
        if fuzzy_equals(max_score, iy): 
            self.iy_matrix.set_pointers(row, col, pointer_row, pointer_col, 'Iy')


    def find_traceback_start(self):
        """
        Finds the location to start the traceback..
        Think carefully about how to set this up for local 

        Returns:
            (max_val, max_loc) where max_val is the best score
            max_loc is a list [] containing tuples with the (i,j) location(s) to start the traceback
             (ex. [(1,2), (3,4)])
        """
        len_a = len(self.align_params.seq_a)
        len_b = len(self.align_params.seq_b)
        starting_pos = (1, 1, '')

        max_loc = set()
        max_value = -np.inf
        if self.align_params.global_alignment:
            #max in last row
            for i in np.arange(1, len_b + 1):
                m = self.m_matrix.get_score(len_a, i)
                ix = self.ix_matrix.get_score(len_a, i)
                iy = self.iy_matrix.get_score(len_a, i)
                curr_max = max(m, ix, iy) 
                tups = (int(len_a), int(i))
                if fuzzy_equals(max_value, curr_max):
                    max_loc.add(tups)
                elif max_value < curr_max:
                    max_loc = set([tups])
                    max_value = curr_max
           
            #max in last col
            for j in np.arange(1, len_a): #not +1 to avoid duplicate of last element
                m = self.m_matrix.get_score(j, len_b)
                ix = self.ix_matrix.get_score(j, len_b)
                iy = self.iy_matrix.get_score(j, len_b)
                curr_max = max(m, ix, iy) 
                tups = (int(j), int(len_b))

                if fuzzy_equals(max_value, curr_max):
                    max_loc.add(tups)
                elif max_value < curr_max:
                    max_loc = set([tups])
                    max_value = curr_max
        else:
            for i in np.arange(1, len_a + 1):
                for j in np.arange(1, len_b + 1):
                    curr_max = self.m_matrix.get_score(i, j)
                    tups = (int(i), int(j))
                    if fuzzy_equals(max_value, curr_max):
                        max_loc.add(tups)
                    elif max_value < curr_max:
                        max_loc = set([tups])
                        max_value = curr_max

        return (float(max_value), max_loc)

    def recursion(self, matrix, row, col, seq_a, seq_b):
        """
        Recursion method to find the traceback path. If not at the end of the 
        sequence (i = 0 or j = 0 for global, i = 0 or j = 0 or pointers = None for global), 
        adds the appropriate letter or gap to each sequence and recursively calls the 
        next pointer(s) to complete the traceback

        Input:
           matrix = current matrix looked at during traceback
           row = the row index to update
           col = the column index to update
           seq_a = current A sequence
           seq_b = current B sequence

        Returns:
            A list of completed traceback sequences
        """
        final_seqs = []

        if row == 0 or col == 0:
            return [(seq_a, seq_b)]

        if matrix == 'M':
            pointers = self.m_matrix.get_pointers(row, col)
            if not pointers and not self.align_params.global_alignment:
                return [(seq_a, seq_b)]
            seq_a = self.align_params.seq_a[row - 1] + seq_a
            seq_b = self.align_params.seq_b[col - 1] + seq_b
        elif matrix == 'Ix':
            pointers = self.ix_matrix.get_pointers(row, col)
            if not pointers and not self.align_params.global_alignment:
                return [(seq_a, seq_b)]
            seq_a = self.align_params.seq_a[row - 1] + seq_a
            seq_b = '_' + seq_b
        else:
            pointers = self.iy_matrix.get_pointers(row, col)
            if not pointers and not self.align_params.global_alignment:
                return [(seq_a, seq_b)]
            seq_a = '_' + seq_a
            seq_b = self.align_params.seq_b[col- 1] + seq_b
        
        for pointer in pointers:
            point_row, point_col, matrix = pointer
            final_seqs += self.recursion(matrix, point_row, point_col, seq_a, seq_b)

        return final_seqs


    def traceback(self):
        """
        Performs a traceback.

        Returns:
            A the rounded alignment score and list of completed traceback sequences
        """
        final_seqs = []
        starting_pos = self.find_traceback_start()
        max_val = starting_pos[0]
        max_positions = starting_pos[1]

        if self.align_params.global_alignment:
            for position in max_positions:
                x_pos = position[0]
                y_pos = position[1]
                for matrix in [self.m_matrix, self.ix_matrix, self.iy_matrix]:
                    if fuzzy_equals(self.m_matrix.get_score(x_pos, y_pos), max_val):
                        final_seqs += self.recursion('M', x_pos, y_pos, "", "")
                    if fuzzy_equals(self.ix_matrix.get_score(x_pos, y_pos), max_val):
                        final_seqs += self.recursion('Ix', x_pos, y_pos, "", "")
                    if fuzzy_equals(self.iy_matrix.get_score(x_pos, y_pos), max_val):
                        final_seqs += self.recursion('Iy', x_pos, y_pos, "", "")
        else:
            for position in max_positions:
                x_pos = position[0]
                y_pos = position[1]
                final_seqs += self.recursion('M', x_pos, y_pos, "", "")

        unique_pairs = list(set(final_seqs))
        return round(float(max_val), 1), unique_pairs

                

    def write_output(self):
        """
        Method to write the output to the output file.

        Output:
            An output file with the alignment sequences written
        """
        f = open(self.output_file, 'w+')
        f.write("{}\n".format(self.align_score))
        f.write("\n")
        for pair in self.alignments:
            f.write("{}\n{}\n".format(pair[0], pair[1]))
            f.write("\n")
        f.close()





def main():

    # check that the file is being properly used
    if (len(sys.argv) !=3):
        print("Please specify an input file and an output file as args.")
        return
        
    # input variables
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # create an align object and run
    align = Align(input_file, output_file)
    align.align()


if __name__=="__main__":
    main()
