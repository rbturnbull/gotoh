# cython: boundscheck=False, wraparound=False
# cython: cdivision=False
# cython: c_string_type=unicode, c_string_encoding=utf8
# cython: language_level=3
import numpy as np
cimport numpy as np

import faulthandler; faulthandler.enable()


cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
cdef extern from "Python.h":
    ctypedef void PyObject
    PyObject *PyString_FromStringAndSize(char *, size_t)
    int _PyString_Resize(PyObject **, size_t)
    char * PyUnicode_AsUTF8(PyObject *)


cimport cython

ctypedef np.int_t DTYPE_t

ctypedef np.int_t DTYPE_INT
ctypedef np.uint_t DTYPE_UINT
ctypedef np.int8_t DTYPE_BOOL

cdef size_t UP = 1, LEFT = 2, DIAG = 3, NONE = 4


cpdef (int, int, int, int) counts(str seqi, str seqj,
                              float match = 0,
                              float mismatch = -1,
                              float gap_open = -1,
                              float gap_extend = -1):
    """
    Perform a global sequence alignment (Needleman-Wunsch) with affine gap penalties (Gotoh) on seqi and seqj. It returns a tuple with the number of characters that match, mismatch, open gaps and extend gaps.
    The scores/penalties are given as arguments and the defaults correspond to basic Levenshtein distance.

    Based almost entirely on Brent Pedersen’s nwalign: https://bitbucket.org/brentp/biostuff/
    """
    cdef bint flip = 0
    
    cdef size_t max_j = len(seqj) 
    cdef size_t max_i = len(seqi) 
    if max_i == max_j == 0:
        return 0, 0, 0, 0


    if max_j > max_i:
        flip = 1
        seqi, seqj = seqj, seqi
        max_i, max_j = max_j, max_i

    # need to initialize j for the case when it's a zero-length string.
    cdef size_t i = 0, j = 0, p
    cdef float diag_score, up_score, left_score, tscore

    cdef char ci, cj
    cdef int zero=0, one=1
    cdef int match_count=0, mismatch_count=0, gap_count=0, extension_count=0

    assert gap_extend <= 0, "gap_extend penalty must be <= 0"
    assert gap_open <= 0, "gap_open must be <= 0"

    cdef np.ndarray[DTYPE_BOOL, ndim=1] agap_i = np.ones((max_i + 1,), dtype=np.int8)
    cdef np.ndarray[DTYPE_BOOL, ndim=1] agap_j = np.ones((max_j + 1,), dtype=np.int8)

    cdef np.ndarray[float, ndim=2] score = np.empty((max_i + 1, max_j + 1), dtype=np.float32)
    cdef np.ndarray[DTYPE_UINT, ndim=2] pointer = np.empty((max_i + 1, max_j + 1), dtype=np.uint)
  
    pointer[0, 0] = NONE
    score[0, 0] = 0
    
    pointer[0, 1:] = LEFT
    pointer[1:, 0] = UP

    score[0, 1:] = gap_open + gap_extend * np.arange(0, max_j, dtype=np.float32)
    score[1:, 0] = gap_open + gap_extend * np.arange(0, max_i, dtype=np.float32)

    agap_i[0] = zero

    for i in range(1, max_i + 1):
        ci = seqi[<size_t>(i - 1)]

        agap_j[0] = zero
        for j in range(1, max_j + 1):
            agap_j[j] = one
            cj = seqj[<size_t>(j - 1)]
            tscore = match if ci == cj else mismatch

            diag_score = score[<size_t>(i - 1), <size_t>(j - 1)] + tscore
            up_score   = score[<size_t>(i - 1), j] + (gap_open if agap_i[<size_t>(i - 1)] == zero else gap_extend)
            left_score = score[i, <size_t>(j - 1)] + (gap_open if agap_j[<size_t>(j - 1)] == zero else gap_extend)

            """
            fix cases where scores are tied.
            choose diagonal when not at the ends of either string.
            and choose up/left (gap) when at the end of the string.
            this forces gaps to the ends of the aligments.
            """
            if diag_score == left_score:
                # so here. if we're at the end we choose a gap.
                # otherwise, we choose diag
                if i == max_i or i == 1:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT
                else: # we want a diagonal
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero
            elif diag_score == up_score:
                if j == max_j or j == 1:
                    score[i, j] = up_score
                    pointer[i, j] = UP
                else: # we want a diagonal
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero
            # end of ambiguous score checks.

            elif up_score > diag_score: #checked.
                if up_score > left_score:
                    score[i, j]  = up_score
                    pointer[i, j] = UP
                else:
                    score[i, j]   = left_score
                    pointer[i, j] = LEFT
            elif diag_score > up_score:
                if left_score > diag_score:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT
                else:
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero

    p = pointer[i, j]
    cdef size_t previous_pointer = NONE
    while p != NONE:
        if p == DIAG:
            i -= 1
            j -= 1
            if seqj[j] == seqi[i]:
                match_count += 1
            else:
                mismatch_count += 1
        elif p == LEFT:
            j -= 1
            if p == previous_pointer:
                extension_count += 1
            else:
                gap_count += 1
            
        elif p == UP:
            i -= 1
            if p == previous_pointer:
                extension_count += 1
            else:
                gap_count += 1
        else:
            raise Exception('Bad Pointer!:pointer: %i', p)
        previous_pointer = p
        p = pointer[i, j]
        
        
    return match_count, mismatch_count, gap_count, extension_count
    

cpdef (float) score(str seqi, str seqj,
                              float match = 0,
                              float mismatch = -1,
                              float gap_open = -1,
                              float gap_extend = -1):
    """
    Returns the score of a sequence alignment.
    """
    my_counts = counts(seqi,seqj,match,mismatch,gap_open,gap_extend)
    return my_counts[0]*match + my_counts[1]*mismatch + my_counts[2]*gap_open + my_counts[3]*gap_extend


cpdef (int) nonmatches(str seqi, str seqj,
                              float match = 0,
                              float mismatch = -1,
                              float gap_open = -1,
                              float gap_extend = -1):
    """
    Returns the number of different characters that do not match.
    """
    my_counts = counts(seqi,seqj,match,mismatch,gap_open,gap_extend)
    return my_counts[1] + my_counts[2] + my_counts[3]


cpdef (float) weighted_nonmatches(str seqi, str seqj,
                              float match = 0,
                              float mismatch = -1,
                              float gap_open = -1,
                              float gap_extend = -1):
    """
    Returns the number of different characters that do not match, multiplied by the negative weight for each category.
    """
    my_counts = counts(seqi,seqj,match,mismatch,gap_open,gap_extend)
    return -(my_counts[1]*mismatch + my_counts[2]*gap_open + my_counts[3]*gap_extend)



@cython.boundscheck(False)
@cython.nonecheck(False)
def msa(np.ndarray[DTYPE_t, ndim=2] seqj, np.ndarray[DTYPE_t, ndim=2] seqi, float gap_open=-1, float gap_extend=-1, object matrix=None, int visualize=0, bint show_score=0):
    """
    perform a global sequence alignment (needleman-wunsch) for multiple sequences
    """
    cdef bint flip = 0

    # if seqj == None:
    #     raise ValueError(f"seqj is null")
    # if seqi == None:
    #     raise ValueError(f"seqi is null")

    cdef size_t max_j = seqj.shape[0]
    cdef size_t max_i = seqi.shape[0]


    if max_i == max_j == 0:
        return np.empty( shape=(0, seqi.shape[1]+seqj.shape[1] ), dtype=int )

    if max_j > max_i:
        flip = 1
        seqi, seqj = seqj, seqi
        max_i, max_j = max_j, max_i


    # need to initialize j for the case when it's a zero-length string.
    cdef size_t i = 0, j = 0, seqlen, p
    cdef float diag_score, up_score, left_score, tscore

    #cdef np.ndarray[DTYPE_t, ndim=1] ci, cj
    cdef int zero=0, one=1

    assert gap_extend <= 0, "gap_extend penalty must be <= 0"
    assert gap_open <= 0, "gap_open must be <= 0"

    cdef np.ndarray[DTYPE_BOOL, ndim=1] agap_i = np.ones((max_i + 1,), dtype=np.int8)
    cdef np.ndarray[DTYPE_BOOL, ndim=1] agap_j = np.ones((max_j + 1,), dtype=np.int8)

    cdef np.ndarray[np.float32_t, ndim=2] score = np.empty((max_i + 1, max_j + 1), dtype=np.float32)
    cdef np.ndarray[DTYPE_UINT, ndim=2] pointer = np.empty((max_i + 1, max_j + 1), dtype=np.uint)

    cdef np.ndarray[np.float32_t, ndim=2] scoring_matrix = matrix
  

    pointer[0, 0] = NONE
    score[0, 0] = 0
    
    pointer[0, 1:] = LEFT
    pointer[1:, 0] = UP

    score[0, 1:] = gap_open + gap_extend * np.arange(0, max_j, dtype=np.float32)
    score[1:, 0] = gap_open + gap_extend * np.arange(0, max_i, dtype=np.float32)

    agap_i[0] = zero
    cdef int GAP = -1

    for i in range(1, max_i + 1):
        ci = set(seqi[<size_t>(i - 1)])

        agap_j[0] = zero
        for j in range(1, max_j + 1):
            agap_j[j] = one
            cj = set(seqj[<size_t>(j - 1)])

            tscore = 0.0
            for token_i in ci:
                for token_j in cj:
                    if token_i == GAP:
                        # If both gaps
                        if token_j == GAP:
                            myscore = 0
                        else:
                            myscore = gap_open if agap_i[<size_t>(i - 1)] == zero else gap_extend
                    elif token_j == GAP:
                        myscore = gap_open if agap_j[<size_t>(j - 1)] == zero else gap_extend
                    elif matrix is None:
                        myscore = 1 if token_i == token_j else -1
                    elif token_i >= token_j:
                        myscore = scoring_matrix[ token_i,token_j ]
                    else:
                        myscore = scoring_matrix[ token_j,token_i ]
                    
                    tscore += myscore
            tscore /= len(ci)*len(cj) # See eq 8 in Spencer and Howe, 2004 (p. 262)


            diag_score = score[<size_t>(i - 1), <size_t>(j - 1)] + tscore
            up_score   = score[<size_t>(i - 1), j] + (gap_open if agap_i[<size_t>(i - 1)] == zero else gap_extend)
            left_score = score[i, <size_t>(j - 1)] + (gap_open if agap_j[<size_t>(j - 1)] == zero else gap_extend)

            """
            fix cases where scores are tied.
            choose diagonal when not at the ends of either string.
            and choose up/left (gap) when at the end of the string.
            this forces gaps to the ends of the aligments.
            """
            if diag_score == left_score:
                # so here. if we're at the end we choose a gap.
                # otherwise, we choose diag
                if i == max_i or i == 1:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT
                else: # we want a diagonal
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero
            elif diag_score == up_score:
                if j == max_j or j == 1:
                    score[i, j] = up_score
                    pointer[i, j] = UP
                else: # we want a diagonal
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero
            # end of ambiguous score checks.

            elif up_score > diag_score: #checked.
                if up_score > left_score:
                    score[i, j]  = up_score
                    pointer[i, j] = UP
                else:
                    score[i, j]   = left_score
                    pointer[i, j] = LEFT
            elif diag_score > up_score:
                if left_score > diag_score:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT
                else:
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero

    seqlen = max_i + max_j
    if visualize:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        im = ax.imshow(score, cmap=plt.get_cmap('hot'))
        fig.colorbar(im)
        plt.show()

    if show_score:
        print("Scores:")
        print(score)

    cdef size_t depth_i =  seqi.shape[1]
    cdef size_t depth_j =  seqj.shape[1]

    cdef int score_max #, = score[:, -1].max()
    cdef np.ndarray[DTYPE_INT, ndim=2] alignment = np.empty((seqlen, depth_i + depth_j), dtype=np.int)
    
    cdef int alignment_index = seqlen - 1
    
    p = pointer[i, j]

    if flip:
        indexes_i = np.arange(depth_i)
        indexes_j = np.arange(depth_i,depth_i+depth_j)
    else:
        indexes_j = np.arange(depth_j)
        indexes_i = np.arange(depth_j,depth_i+depth_j)

    while p != NONE:
        if p == DIAG:
            i -= 1
            j -= 1

            alignment[alignment_index,indexes_i] = seqi[i,:]
            alignment[alignment_index,indexes_j] = seqj[j,:]
        elif p == LEFT:
            j -= 1

            alignment[alignment_index,indexes_i] = GAP
            alignment[alignment_index,indexes_j] = seqj[j,:]
        elif p == UP:
            i -= 1
            alignment[alignment_index,indexes_i] = seqi[i,:]
            alignment[alignment_index,indexes_j] = GAP
        else:
            raise Exception('Error with pointer: %i', p)
        alignment_index -= 1
        p = pointer[i, j]

    return alignment[alignment_index+1:,:]



@cython.boundscheck(False)
@cython.nonecheck(False)
def pointers(np.ndarray[DTYPE_t, ndim=2] seqj, np.ndarray[DTYPE_t, ndim=2] seqi, float gap_open=-1, float gap_extend=-1, object matrix=None, int visualize=0, bint show_score=0):
    """
    perform a global sequence alignment (needleman-wunsch) for multiple sequences.
    Returns an array with the pointers to generate the alignment.
    """
    cdef bint flip = 0

    # if seqj == None:
    #     raise ValueError(f"seqj is null")
    # if seqi == None:
    #     raise ValueError(f"seqi is null")

    cdef size_t max_j = seqj.shape[0]
    cdef size_t max_i = seqi.shape[0]


    if max_i == max_j == 0:
        return np.empty( shape=(0, seqi.shape[1]+seqj.shape[1] ), dtype=int )

    if max_j > max_i:
        flip = 1
        seqi, seqj = seqj, seqi
        max_i, max_j = max_j, max_i


    # need to initialize j for the case when it's a zero-length string.
    cdef size_t i = 0, j = 0, seqlen, p
    cdef float diag_score, up_score, left_score, tscore

    #cdef np.ndarray[DTYPE_t, ndim=1] ci, cj
    cdef int zero=0, one=1

    assert gap_extend <= 0, "gap_extend penalty must be <= 0"
    assert gap_open <= 0, "gap_open must be <= 0"

    cdef np.ndarray[DTYPE_BOOL, ndim=1] agap_i = np.ones((max_i + 1,), dtype=np.int8)
    cdef np.ndarray[DTYPE_BOOL, ndim=1] agap_j = np.ones((max_j + 1,), dtype=np.int8)

    cdef np.ndarray[np.float32_t, ndim=2] score = np.empty((max_i + 1, max_j + 1), dtype=np.float32)
    cdef np.ndarray[DTYPE_UINT, ndim=2] pointer = np.empty((max_i + 1, max_j + 1), dtype=np.uint)

    cdef np.ndarray[np.float32_t, ndim=2] scoring_matrix = matrix
    cdef int GAP = -1  

    pointer[0, 0] = NONE
    score[0, 0] = 0
    
    pointer[0, 1:] = LEFT
    pointer[1:, 0] = UP

    score[0, 1:] = gap_open + gap_extend * np.arange(0, max_j, dtype=np.float32)
    score[1:, 0] = gap_open + gap_extend * np.arange(0, max_i, dtype=np.float32)

    agap_i[0] = zero

    for i in range(1, max_i + 1):
        ci = set(seqi[<size_t>(i - 1)])

        agap_j[0] = zero
        for j in range(1, max_j + 1):
            agap_j[j] = one
            cj = set(seqj[<size_t>(j - 1)])

            tscore = 0.0
            for token_i in ci:
                for token_j in cj:
                    if token_i == GAP:
                        # If both gaps
                        if token_j == GAP:
                            myscore = 0
                        else:
                            myscore = gap_open if agap_i[<size_t>(i - 1)] == zero else gap_extend
                    elif token_j == GAP:
                        myscore = gap_open if agap_j[<size_t>(j - 1)] == zero else gap_extend
                    elif matrix is None:
                        myscore = 1 if token_i == token_j else -1
                    elif token_i >= token_j:
                        myscore = scoring_matrix[ token_i,token_j ]
                    else:
                        myscore = scoring_matrix[ token_j,token_i ]
                    
                    tscore += myscore
            tscore /= len(ci)*len(cj) # See eq 8 in Spencer and Howe, 2004 (p. 262)


            diag_score = score[<size_t>(i - 1), <size_t>(j - 1)] + tscore
            up_score   = score[<size_t>(i - 1), j] + (gap_open if agap_i[<size_t>(i - 1)] == zero else gap_extend)
            left_score = score[i, <size_t>(j - 1)] + (gap_open if agap_j[<size_t>(j - 1)] == zero else gap_extend)

            """
            fix cases where scores are tied.
            choose diagonal when not at the ends of either string.
            and choose up/left (gap) when at the end of the string.
            this forces gaps to the ends of the aligments.
            """
            if diag_score == left_score:
                # so here. if we're at the end we choose a gap.
                # otherwise, we choose diag
                if i == max_i or i == 1:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT
                else: # we want a diagonal
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero
            elif diag_score == up_score:
                if j == max_j or j == 1:
                    score[i, j] = up_score
                    pointer[i, j] = UP
                else: # we want a diagonal
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero
            # end of ambiguous score checks.

            elif up_score > diag_score: #checked.
                if up_score > left_score:
                    score[i, j]  = up_score
                    pointer[i, j] = UP
                else:
                    score[i, j]   = left_score
                    pointer[i, j] = LEFT
            elif diag_score > up_score:
                if left_score > diag_score:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT
                else:
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero

    if show_score:
        print("----------------------")
        print("Scores:")
        print(score)
        print("^^^^^^^^^^^^^^^^^^^^^^")

    return pointer


def pointer_constants():
    return UP, LEFT, DIAG, NONE

def pointers_ascii(pointers):
    pointer_to_char = {
        UP: "^",
        LEFT: "<",
        DIAG: "\\",
        NONE: ".",
    }
    rows, columns = pointers.shape
    string = ""
    for i in range(rows):
        for j in range(columns):
            pointer = pointers[i,j]
            string += " " + pointer_to_char[ pointer ] + " "


        string += "\n"
    return string