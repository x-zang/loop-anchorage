"""
Author: Carl Zang
Date Started: Jul 2023
Last updated: Aug 17 2023 or check git commit 

Version:    v0.3.0

Description:
    The main pipeline of anchorage. It performs contig assembly for LoopSeq SLR.

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

from anchorGuidedAssembler import AnchorGuidedAssembler
from contigStatistician import ContigStatistician
from aligner import PseudoAligner
import aligner
from loop.common.utils import RunShellCommand

from sys import argv
import sys
import argparse
import os

def anchorage(args):
    execute_spades(args)
    spades_gfa = os.path.join(args.output_prefix, "assembly_graph_with_scaffolds.gfa")
    execute_AGA(args, spades_gfa)
    error_correction()
    print("Anchorage completed running!")
    return 0

def execute_spades(args):
    cmd_spades = []
    cmd_spades.append("spades.py")
    cmd_spades.append("--pe-1  0 {}".format(args.read1_fq))
    cmd_spades.append("--pe-2  0 {}".format(args.read2_fq))
    cmd_spades.append("-k  {}".format(args.k))
    cmd_spades.append("-o {}".format(args.output_prefix))
    cmd_spades.append("--careful --sc -t 8")       
    cmd_spades.append("--phred-offset  33")
    cmd_spades.append("--disable-gzip-output")
    cmd_line = " ".join(cmd_spades)
    RunShellCommand(cmd_line, 'Calling spades')
    return 0

def execute_AGA(args, gfa):
    """
        output file: args.output_prefix + ".fa", default "anchorage_contig.fa"
    """
    AGA = AnchorGuidedAssembler(
                        # args.gfa,
                        gfa,
                        args.anchor_start,
                        args.anchor_end,
                        left_fq=args.read1_fq,
                        right_fq=args.read2_fq,
                        output_contig_file = args.output_prefix,
                        ht_index_1=args.ht_index_1,
                        ht_index_2=args.ht_index_2,
                        barcode_trim = args.trim_barcode,
                        max_nm_anchors = args.max_nm_anchors,
                        barcode_len = args.contig_barcode_len,   
                        verbose=args.verbose)
    return 0

def error_correction():
    #TODO: re-align reads back to contig and perform error corrections
    pass

def execute_argparse():
    parser = argparse.ArgumentParser()
    requiredarg = parser.add_argument_group('Required arguments')
    recommendarg = parser.add_argument_group('Recommended arguments')
    configarg = parser.add_argument_group('Algorithm configuration arguments')

    # requiredarg.add_argument('-i', '--gfa', required=True, metavar='\b',
    #                             help="GFA file produced by various modes of SPAdes.")
    requiredarg.add_argument('-s1', '--anchor_start', required=True, metavar='\b')
    requiredarg.add_argument('-s2', '--anchor_end', required=True, metavar='\b')
    requiredarg.add_argument('-r1', '--read1_fq', required=True, metavar='\b')
    requiredarg.add_argument('-r2', '--read2_fq', required=True, metavar='\b')

    recommendarg.add_argument('-h1', '--ht_index_1', default="", metavar='\b',
                        help="A variable start anchor following anchor_start immediately in molecule, e.g. HT_index_1 or equivalent")
    recommendarg.add_argument('-h2', '--ht_index_2', default="", metavar='\b',
                        help="A variable end anchor following anchor_end immediately in molecule, e.g. HT_index_2 or equivalent")   
    recommendarg.add_argument('-o', '--output_prefix', default="anchorage_contig", required=False, metavar='\b',
                        help='output file prefix, default: anchor_guide_contig')
    recommendarg.add_argument('-k', default="21,33,55,77,99", required=False, metavar='\b',
                        help='a series of k-mer sizes, default 21,33,55,77,99.')
    
    configarg.add_argument('--contig_barcode_len', type=int, metavar='\b',
                        help='length of contig barcode')
    
    configarg.add_argument('--trim_barcode', type=bool, default=True, metavar='\b',
                        help='if trimming barcode from contig ')

    configarg.add_argument('--max_nm_anchors', type=int, default=2, metavar='\b',
                        help='maximum number of anchors permitted')

    configarg.add_argument('--use_ambiguous_anchor', type=bool, metavar='\b',
                        help='if anchor node is ambiguous, use the max-depth one')
    # configarg.add_argument('--allow_repeats_same_orientation', type=bool, metavar='\b')
    # configarg.add_argument('--allow_repeats_revcomp', type=bool, metavar='\b')
    # configarg.add_argument('--always_keep_ndoe_longer_than_length_good', type=bool, metavar='\b')
    configarg.add_argument('--verbose', type=bool, metavar='\b')
    
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = execute_argparse()
    anchorage(args)
