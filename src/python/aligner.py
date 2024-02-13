"""
Author: Carl Zang
Date Started: Jul 2023
Last updated: Jul 25 2023 or check git history 

Version:    v0.0.1


BSD 3-Clause License

Copyright (c) 2024, Element Biosciences 

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


"""

# external library
import numpy as np
# standard library
import sys
import os
import math



def smith_waterman(seq1, seq2, match_score=1, mismatch_penalty=-1, indel_penalty=-1):
    if not set(seq1 + seq2).issubset(set("ATCG")):
        print("WARNING:\t RealAligner.smith_waterman: sequenes have non-ATCG bases. Proceed anyway.", file=sys.stdout)

    mx = np.zeros(((len(seq1) + 1), (len(seq2) + 1)), dtype=int)
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            examine_match = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty
            consume_both = mx[i-1][j-1] + examine_match
            consume_seq1 = mx[i  ][j-1] + indel_penalty  
            consume_seq2 = mx[i-1][j  ] + indel_penalty
            mx[i][j] = max(0, consume_both, consume_seq1, consume_seq2)
    return max([max(i) for i in mx])


def needleman_wunsch(seq1, seq2, match_score = 1, mismatch_penalty = -1, indel_penalty = -1):
    if not set(seq1 + seq2).issubset(set("ATCG")):
        print("WARNING:\t RealAligner.smith_waterman: sequenes have non-ATCG bases. Proceed anyway.", file=sys.stdout)

    mx = np.zeros(((len(seq1) + 1), (len(seq2) + 1)), dtype=int)
    for i in range(1, len(seq1) + 1):
        mx[i][0] = indel_penalty * i
    for i in range(1, len(seq2) + 1):
        mx[0][i] = indel_penalty * i

    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            examine_match = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty
            consume_both = mx[i-1][j-1] + examine_match
            consume_seq1 = mx[i  ][j-1] + indel_penalty  
            consume_seq2 = mx[i-1][j  ] + indel_penalty
            mx[i][j] = max(consume_both, consume_seq1, consume_seq2)
    return mx[len(seq1)][len(seq2)]


def edit_distance_global(seq1, seq2):
    """
        return number of operations (substitution, insertion, deletion) needed to
            transform seq1 to seq2
    """
    return -needleman_wunsch(seq1, seq2, match_score=0, mismatch_penalty=-1, indel_penalty=-1)
    

def is_subseq(seq1, seq2, nm=0):
    """
        return bool if using number of operations (substitution, insertion, deletion) can
            transform seq1 to a subsequence of seq2
        Assume seq1 is shorter seq2
    """
    if not set(seq1 + seq2).issubset(set("ATCG")):
        print("WARNING:\t sequenes have non-ATCG bases. Proceed anyway.", file=sys.stdout)

    if len(seq1) > len(seq2):
        return False
    
    if nm == 0:
        return (seq1 in seq2)

    mx = np.zeros(((len(seq1) + 1), (len(seq2) + 1)), dtype=int)
    
    indel_penalty = 1
    for i in range(1, len(seq1) + 1):
        mx[i][0] = i * indel_penalty
    # for i in range(1, len(seq2) + 1):
    #     mx[0][i] = i * indel_penalty
    
    for j in range(1, len(seq2) + 1):
        for i in range(1, len(seq1) + 1):
            examine_match = 0 if seq1[i - 1] == seq2[j - 1] else 1
            consume_both = mx[i-1][j-1] + examine_match
            consume_seq1 = mx[i  ][j-1] + 1  
            consume_seq2 = mx[i-1][j  ] + 1 
            mx[i][j] = min(consume_both, consume_seq1, consume_seq2)
        if mx[i][j] <= nm:
            return True
    return False


def subseq_pos(seq1, seq2, nm=0):
    """
        return position of seq2's first subseq (return -1 if unfeasible) 
        if using number of operations (substitution, insertion, deletion) can
            transform seq1 to this subsequence of seq2
        Assume seq1 is shorter seq2
    """
    if not set(seq1 + seq2).issubset(set("ATCG")):
        print("WARNING:\t sequenes have non-ATCG bases. Proceed anyway.", file=sys.stdout)

    if len(seq1) > len(seq2):
        return -1
    
    if nm == 0:
        for i in range(len(seq2) - len(seq1) + 1):
            if seq2[i:].startswith(seq1):
                return i + len(seq1)
        return -1

    mx = np.zeros(((len(seq1) + 1), (len(seq2) + 1)), dtype=int)

    indel_penalty = 1
    for i in range(1, len(seq1) + 1):
        mx[i][0] = i * indel_penalty
    # for i in range(1, len(seq2) + 1):
    #     mx[0][i] = i * indel_penalty

    for j in range(1, len(seq2) + 1):
        for i in range(1, len(seq1) + 1):
            examine_match = 0 if seq1[i - 1] == seq2[j - 1] else 1
            consume_both = mx[i-1][j-1] + examine_match
            consume_seq1 = mx[i  ][j-1] + 1  
            consume_seq2 = mx[i-1][j  ] + 1 
            mx[i][j] = min(consume_both, consume_seq1, consume_seq2)
        if mx[i][j] <= nm:
            return j
    return -1


class PseudoAligner():
    """
        Input:
            Sequence files such as GFA, fasta

        Description:
            Performs pseudo-alignment with kallisto
            May perform exact alignment       
    """
    def __init__():
        raise NotImplementedError
    
    def index(self):
        pass
    
    def get_reads_with_subseq(self):
        pass

    def pseudo_align_reads(self):
        pass

    def pseudo_align_per_read(self):
        pass

    def exact_align_reads(self):
        pass

    def exact_align_per_read(self):
        # run SW, return alignment or CIGAR
        pass
    
    def align_count(self):
        pass


def __debug_tests():
    print("Running tests for PsudoAligner and RealAligner")
    assert smith_waterman('AAAAAA', 'TTTTTT') == 0
    assert smith_waterman('AAAAAA', 'A') == 1
    assert smith_waterman('A', 'AAAAAA') == 1
    assert smith_waterman('ATAACG', 'GCAATA') == 3
    assert smith_waterman('AAATTT', 'TTTCCC', 2, -1, -1) == 6
    assert smith_waterman('ATATAA', 'TTTTAT', 2, 0, 0) == 6
    assert smith_waterman('CGCGCG', 'CGACTCG', 5, -1, -1) == 25 - 2
    assert needleman_wunsch('AAAAAA', 'TTTTTT') == -6
    assert needleman_wunsch('AAAAAA', 'A') == -4
    assert needleman_wunsch('A', 'AAAAAA') == -4
    assert needleman_wunsch('ATAACG', 'GCAATA') == -2
    assert needleman_wunsch('AAATTT', 'TTTCCC', 2, -1, -1) == 0
    assert needleman_wunsch('ATATAA', 'TTTTAT', 2, 0, 0) == 6
    assert needleman_wunsch('CGCGCG', 'CGACTCG', 5, -1, -1) == 25 - 2

    assert edit_distance_global('AAAAAA', 'TTTTTT') == 6
    assert edit_distance_global('AAAAAA', 'A') == 5
    assert edit_distance_global('A', 'AAAAAA') == 5
    assert edit_distance_global('ATAACG', 'GCAATA') == 4
    assert edit_distance_global('AAATTT', 'TTTCCC') == 6
    assert edit_distance_global('ATATAA', 'TTTTAT') == 3
    assert edit_distance_global('CGCGCG', 'CGACTCG') == 2

    assert not is_subseq('AAAAAA', 'TTTTTT') 
    assert not is_subseq('AAAAAA', 'A') 
    # assert is_subseq('A', 'AAAAAA') 
    assert not is_subseq('ATAACG', 'GCAATA') 
    assert is_subseq('TTT', 'AGGGATTTCCC') 
    assert is_subseq('ATATAA', 'ATATAATTTTAT') 
    assert is_subseq('ATATAA', 'TTTTATATATAA')
    assert not is_subseq('CGCGCG', 'CGACTCG', 0)
    assert not is_subseq('CGCGCG', 'CGACTCG', 1)
    assert is_subseq('CGCGCG', 'CGACTCG', 2)
    assert not is_subseq('CGCGCG', 'CGACTCG')
    assert is_subseq('CGAACT', 'CGACTCG', 1)
    assert is_subseq('CGCT', 'CGACTCG', 1)
    assert is_subseq('AGCT', 'CGACTCG', 2)


    assert subseq_pos('AAAAAA', 'TTTTTT') == -1 
    assert subseq_pos('AAAAAA', 'A') == -1
    assert subseq_pos('A', 'AAAAAA') == 1
    assert subseq_pos('C', 'AACAAA') == 3
    assert subseq_pos('ATAACG', 'GCAATA') == -1
    assert subseq_pos('TTT', 'AGGGATTTCCC') == 8
    assert subseq_pos('ATATAA', 'ATATAATTTTAT') == 6
    assert subseq_pos('ATATAA', 'TTTTATATATAA') == 12
    assert subseq_pos('CGCGCG', 'CGACTCG', 2) == 7
    assert subseq_pos('CGAACT', 'CGACTCG', 1) == 5
    assert subseq_pos('CGAACT', 'CGACTCG', 0) == -1
    assert subseq_pos('CGCT', 'CGACTCG', 1) == 5
    assert subseq_pos('CGCT', 'CGACTCGCGCT', 1) == 5
    assert subseq_pos('AAAAAA', 'AAAAAA', 1) == 5
    assert subseq_pos('ATCGTG', 'ATCGTG', 1) == 5
    assert subseq_pos('ATCGTG', 'ATCGTG', 0) == 6
    assert subseq_pos('AGCT', 'CGACTCG', 1) == 5
    assert subseq_pos('AAAAAA', 'CCAAAAAA', 1) == 7
    assert subseq_pos('ATCGTG', 'CCATCGTG', 1) == 7
    assert subseq_pos('ATCGTG', 'CCATCGTG', 0) == 8
    assert subseq_pos('AGCT', 'CCCGACTCG', 1) == 7
    print("Aligner all tests passed.")
    return 0


if __name__ == "__main__":
    if sys.argv[1] in ['--test', '-test', 'test'] and len(sys.argv) == 2:
        __debug_tests()
    else:
        raise NotImplementedError("You should not run PseudoAligner from main. Import and call the object.!")