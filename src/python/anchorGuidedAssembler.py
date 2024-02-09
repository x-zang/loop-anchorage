"""
Author: Carl Zang
Date Started: Jul 2023
Last updated: Aug 17 2023 or check git commit 

Version:    v0.3.0

Help information: run `python anchorGuidedAssembler.py -h`


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
import networkx as nx

# internal and standard library
from contigStatistician import ContigStatistician
from loop.analysis.utility import revcomp
import aligner  
import sys
from sys import argv
import argparse
import logging
import json


class AnchorGuidedAssembler():
    """
        Input:
            gfa from spades (Gfa v1.0)

        Description:
            Extract a single walk from spades gfa1 graph, by:
            (1) Utilizing anchor sequences to break loops and guide start/end of walk
                incld. (1.1) remove undefined seq upstream (resp. down) of anchor_start (resp. _end),
                    (1.2) Resolve inverted loop s.t. final contig contains forward anchors only
                            Inverted loop is the major cause where we have truncated contigs
                            Most other algorithms cannot or ad hoc perform operations on inverted loop
                            without being aware of additional information, such as anchors
            (2) Utilizing expected contig length
            (3) Aggressively removing undesired edges by using coverage estimation
    """


    def __init__(
            self, 
            gfafile, 
            anchor_start, 
            anchor_end, 
            left_fq,
            right_fq,
            output_contig_file = "", 
            ht_index_1 = "",
            ht_index_2 = "",
            max_nm_anchors = 2,
            barcode_trim = True,
            barcode_len = 0,
            asssembly_algo_config = dict(), 
            verbose=False):
        
        # config
        self.logger = logging.getLogger('Anchor-Guided Assembly') 
        self.gfa = gfafile
        self.verbose = verbose
        self.left_fq = left_fq
        self.right_fq = right_fq
        self.anchor_start = anchor_start
        self.anchor_end = anchor_end
        self.ht_index_1 = ht_index_1
        self.ht_index_2 = ht_index_2
        self.max_nm_anchors = max_nm_anchors
        self.nm_anchor_start = -1        
        self.nm_anchor_end = -1
        self.outfile = output_contig_file if output_contig_file else "anchor_guide_contig"
        self.assembler_config = {}
        self.barcode_trim = barcode_trim
        if self.barcode_trim:
            assert barcode_len > 0
            self.barcode_len = barcode_len
        self.get_config(asssembly_algo_config)

        # Stat
        self.contig_stats = None
        self.contig_length_max = -1
        self.contig_length_min = -1
        self.contig_length_fav = -1
        self.cov_max_estimate = -1.0
        self.cov_min_estimate = -1.0
        self.cov_fav_estimate = -1.0
        self.get_contig_stats()


        # internal variables
        self.graph = self.gfa_to_nextworkx(self.gfa, self.logger, self.verbose)
        self.anchor_start_nodes = []
        self.anchor_end_nodes = []
        self.source = None
        self.target = None
        self.contig_path = []
        self.contig_seq = ""

        self.assemble()        
        if not self.contig_seq and not self.contig_path:
            raise nx.NetworkXNoPath('No feasible contig is recovered.')
        else:
            self.write_contig()
            self.logger.info('Anchor-guided assembly is completed. Recovered contig in {}'.format(self.outfile + ".fa"))
        
        return None

    def get_config(self, config):
        default_false_config_names =  [
            'allow_repeats_same_orientation',
            'allow_repeats_revcomp', 
            'always_keep_ndoe_longer_than_length_good'
        ]
        default_true_config_names = [
            'use_ambiguous_anchor'
        ]

        for i in config.keys():
            assert i in default_false_config_names or default_true_config_names

        for i in default_false_config_names:
            if not config or not i in config:
                self.assembler_config[i] = False
            else:
                self.assembler_config[i] = config[i]
        for i in default_true_config_names:
            if not config or not i in config:
                self.assembler_config[i] = True
            else:
                self.assembler_config[i] = config[i]

        assert len(self.assembler_config) == len(default_false_config_names) + len(default_true_config_names)
        return 0

    def get_contig_stats(self):
        """
            get upper/lower/estimated contig length and node depth
            Note: covs are no longer w.r.t. lengths since they have different tolerance
                  covs estimations assume strong uniformity in data
        """
        self.contig_stats = ContigStatistician(self.left_fq, self.right_fq, nameprefix = self.outfile + ".AGA.CStat")
        cs = self.contig_stats
        self.contig_length_fav = cs.length_good
        self.contig_length_min = cs.length_lower_bound
        self.contig_length_max = cs.length_upper_bound        
        self.cov_fav_estimate = cs.cov_estimation
        self.cov_min_estimate = max(cs.cov_lower_bound, 0.1 * self.cov_fav_estimate)
        self.cov_max_estimate = min(cs.cov_upper_bound, 5   * self.cov_fav_estimate)
        msg = "Contig statistics:\n"
        msg += "contig_length_max = {}\n".format(self.contig_length_max)        
        msg += "contig_length_min = {}\n".format(self.contig_length_min)
        msg += "contig_length_fav = {}\n".format(self.contig_length_fav)
        msg += "cov_max_estimate = {}\n".format(self.cov_max_estimate)
        msg += "cov_min_estimate = {}\n".format(self.cov_min_estimate)
        msg += "cov_fav_estimate = {}\n".format(self.cov_fav_estimate)
        self.logger.debug(msg)
        assert self.contig_length_max >= self.contig_length_fav and self.contig_length_fav >= self.contig_length_min      
        assert self.cov_max_estimate >= self.cov_fav_estimate and self.cov_fav_estimate >= self.cov_min_estimate
        return 0

    @staticmethod
    def gfa_to_nextworkx(filename, logger=None, verbose=False):
        """  
        Args:
            filename (str): GFA1 file from spades

        Returns:
            graph:  nx.Graph() representation of GFA
                    Each S is represented by two nodes, '+' and '-'
                    Each L is represented by two edges, both forward and reversed
        """
        graph = nx.DiGraph()
        with open(filename, 'r') as f:
            lines = [x.strip().split('\t') for x in f.readlines()]
            comment_lines = [x for x in lines if x[0].startswith('#')]
            H_lines = [x for x in lines if x[0] == 'H']
            S_lines = [x for x in lines if x[0] == 'S']
            L_lines = [x for x in lines if x[0] == 'L']
            P_lines = [x for x in lines if x[0] == 'P']
        
        # construct nodes from Segment lines
        # each S is represented in TWO nodes, '+' and '-'
        for sline in S_lines:
            recordtype, name, seq = sline[0:3]
            length = len(seq)
            opt_fields = sline[3:] if len(sline) >=4 else None
            opt = AnchorGuidedAssembler.gfa_optional_field_parser(opt_fields)
            depth = opt['DP']
            graph.add_node(name + '+', length=length, depth=depth, 
                           seq=seq, opt=opt)
            graph.add_node(name + '-', length=length, depth=depth, 
                           seq=revcomp(seq), opt=opt)   
        
        # construct edges from Link lines
        # edges are "stranded"
        for lline in L_lines:
            recordtype, source, sourceOrient, target, targetOrient, overlap = lline[0:6]
            opt_fields = lline[6:] if len(sline) >= 7 else None
            opt = AnchorGuidedAssembler.gfa_optional_field_parser(opt_fields)        
            rev_orientation = {'+': '-', '-': '+'}
            # get edge seq from overlap CIGAR
            # seq_overlap is stranded and matches edge orientation (not source or target orientation)
            assert sourceOrient in {'+', '-'} and targetOrient in {'+', '-'}
            assert overlap[-1] == 'M'  # See catch below about parsing spades GFA
            try:
                len_overlap = int(overlap[:-1])
                seq_overlap = graph.nodes[source + sourceOrient]['seq'][-len_overlap:]
            except TypeError:
                raise TypeError( \
                    "The GFA parser is specific to spades output at this stage." + \
                    "L line overlap is assumed to have \'M\' only, but got {}.", overlap)                 
            if graph.has_edge(source + sourceOrient, target + targetOrient):
                pass 
            else:
                graph.add_edge(source + sourceOrient, target + targetOrient,
                               len_over=len_overlap, seq=seq_overlap, opt=opt)
            if graph.has_edge(target + rev_orientation[targetOrient], source + rev_orientation[sourceOrient]):
                pass
            else:
                graph.add_edge(target + rev_orientation[targetOrient], source + rev_orientation[sourceOrient], 
                            len_over=len_overlap, seq=revcomp(seq_overlap), opt=opt) 
        
        if verbose:
            print("Finished loading GFA to nx.DiGraph, from file {}".format(filename))
            for i in graph.nodes:
                print("Node:",i, graph.nodes[i])
            for i in graph.edges:
                print("Edge:", i, graph.edges[i])

        return graph


    @staticmethod
    def gfa_optional_field_parser(optional_fields):
        """
            input: a list of optional field strings [TAG:TYPE:VALUE]
            Output: a dictionary {TAG:VALUE}
        """
        if not optional_fields:
            return {}
        
        d = {}
        for opt in optional_fields:
            tag, _type, value = opt.strip().split(':', 2)
            assert tag not in d
            if _type == 'A':
                value = str(value)
                assert len(value) <= 1
            elif _type == 'i':
                value = int(value)
            elif _type == 'f':
                value = float(value) 
            elif _type == 'Z':
                value = str(value)
            elif _type == 'J':
                value = json.load(value)
            elif _type == 'H':
                raise Exception('optional field parser did not implement Hex yet')
            elif _type == 'B':
                raise Exception('optional field parser did not implement Array of int/float yet')
            else:
                raise KeyError('optional field has undefined TYPE')
            d[tag] = value
        return d

    def anchor_nodes(self, s, h, nm):
        """
            Input:
                s: anchor sequence
                h: ht_index, assumed to follow anchor immediately
                nm: number of mismatches permitted in total (s+h)
            Returns: Node list containg desired sequence s + h
            Exception: if two nodes contain sequence s
        """
        anchor_nodes = []
        for node in self.graph.nodes():
            node_seq = self.graph.nodes[node]['seq']
            if aligner.is_subseq(s, node_seq, nm) and aligner.is_subseq(s + h, node_seq, nm):
                anchor_nodes.append((node))
                   
        if len(anchor_nodes) == 1:
            return anchor_nodes        
        elif self.assembler_config['use_ambiguous_anchor']:
            return anchor_nodes
        else:
            raise nx.AmbiguousSolution('Anchor {} found in more than one nodes: {}'.\
                                    format(s, anchor_nodes))    

    def anchor_edges(self, s, h):
        """
            Returns: One or none size-2 list [Source Node, Target Node] of an Edge containg desired sequence s
            Exception: if two edges contain sequence s
        """
        raise NotImplementedError("Assuming index length is always smaller than overlap length. No chance of using anchor_edge()")

    def split_anchor_nodes(self):
        """
            Description:
                Breaks graph at start_node_list and end_node_list
                Splitting last (resp. first) node in start_node_list (resp. end) to two nodes, 
                one with 0 in-degree, the other with 0 out-degree
                Note:
                1. Do NOT split their revcomp node
                2. Self-loop will be properly handeled as splitting happens sequentially. The connection is kept.
                    e.g. A- to A+ will be newAhead to A-, newAtail to A+, newAtail to newAhead
        """
        # instead of splitting, make two new nodes instead
        # needs delicate handeling when two nodes forms circulation
        # split start node
        start = self.source
        if not start.startswith('AGAsplitstart_'):
            split_start = 'AGAsplitstart_' + start
            seq = self.graph.nodes[start]['seq']  
            dp = self.graph.nodes[start]['depth']
            pos = aligner.subseq_pos(self.anchor_start, seq, self.nm_anchor_start)
            assert self.nm_anchor_start >= 0
            assert pos >= 0 and pos <= len(seq)
            split_seq = seq[pos:]
            length = len(split_seq)
            self.graph.add_node(split_start, length=length, depth=dp, seq=split_seq)   
            
            for e in self.graph.out_edges(start):
                edge_target = e[1]
                len_over = self.graph.edges[e]['len_over']
                seq_overlap = self.graph.edges[e]['seq']
                # anchor is in overlapping region
                # only overlaping seq will be recorded 
                if len_over >= length: 
                    seq_overlap = split_seq
                self.graph.add_edge(split_start, edge_target, len_over=len_over, seq=seq_overlap) 
            self.source = split_start
            print("Split Node:", split_start, self.graph.nodes[split_start])

        # split end node 
        end = self.target
        if not end.startswith('AGAsplitend_'):
            split_end = 'AGAsplitend_' + end  
            seq = self.graph.nodes[end]['seq']  
            dp = self.graph.nodes[end]['depth']
            pos = len(seq) - aligner.subseq_pos(revcomp(self.anchor_end), revcomp(seq), self.nm_anchor_end)
            assert self.nm_anchor_end >= 0
            assert pos >= 0 
            assert pos <= len(seq)
            split_seq = seq[:pos]
            length = len(split_seq)
            self.graph.add_node(split_end, length=length, depth=dp, seq=split_seq) 
            
            for e in self.graph.in_edges(end):
                edge_source = e[0]
                len_over = self.graph.edges[e]['len_over']
                seq_overlap = self.graph.edges[e]['seq']
                # anchor is in overlapping region
                # only overlaping seq will be recorded 
                if len_over >= length: 
                    seq_overlap = split_seq
                self.graph.add_edge(edge_source, split_end, len_over=len_over, seq=seq_overlap) 
            self.target = split_end
            print("Split Node:", split_end, self.graph.nodes[split_end])

        return 0

    def get_anchor_node_lists(self, mismatch):
        """
            Description:
                Retrieve start anchor nodes list and end anchor nodes list
                At least one s-t path is present
                Attempts the following: (0) whole sequence of anchor + ht_index_1 will be considered
                (1) anchor w/o mismatches (2) anchor with fewer mismatches
        """
        fix_start = self.anchor_start
        fix_end = self.anchor_end
        flex_start = self.ht_index_1
        flex_end = self.ht_index_2
               
        anchor_start_nodes = []
        anchor_end_nodes   = []
        for i in range(mismatch + 1):
            # if len(anchor_start_nodes) == 0:
            self.logger.info('Identifying anchor start nodes permitting {} mismatches'.format(i))
            anchor_start_nodes = self.anchor_nodes(fix_start, flex_start, nm=i)
            self.nm_anchor_start = i

            # if len(anchor_end_nodes) == 0:
            self.logger.info('Identifying anchor end nodes permitting {} mismatches'.format(i))
            anchor_end_nodes = self.anchor_nodes(fix_end, flex_end, nm=i)
            self.nm_anchor_end = i

        self.logger.info('Anchor start nodes: {}'.format(anchor_start_nodes))
        self.logger.info('Anchor end nodes: {}'.format(anchor_start_nodes))
        
        if len(anchor_start_nodes) == 0 or len(anchor_end_nodes) == 0:
            raise nx.NodeNotFound("No anchor node found even if with max {} mismatches".format(mismatch))        

        self.anchor_start_nodes = anchor_start_nodes
        self.anchor_end_nodes = anchor_end_nodes

        return 0
            
    def source_must_reach_target(self):
        assert self.source
        assert self.target
        assert self.graph
        if not nx.has_path(self.graph, self.source, self.target):
            raise nx.NetworkXUnfeasible(
                'No path exists between source and target. Contigs does not exist in graph!')


    @staticmethod
    def node_path_to_seq(graph, node_list):
        """
            Input: 
                a list of names of a consecutive node path
            Return:
                a string of sequence of the path, after removing overlapped bases between nodes   
        """
        if len(node_list) == 0:
            return ''
        if len(node_list) == 1:
            return graph.nodes[node_list[0]]['seq']
        
        # first & last node need to be handeled separately
        seq = ""
        seq += graph.nodes[node_list[0]]['seq']
        for i in range(1, len(node_list)):
            s = node_list[i - 1]
            t = node_list[i]
            try:
                overlap_len = graph.edges[(s,t)]['len_over']
                # overlap_seq = graph.edges[(s,t)]['seq']
            except KeyError:
                raise nx.NodeNotFound('Path {} is not valid or edge {} misses length property'.format(node_list, (s,t)))
            tseq = graph.nodes[t]['seq']
            if len(tseq) >= overlap_len:
                seq += tseq[overlap_len: ]
            else:
                seq = seq[:len(tseq) - overlap_len]
        # last_seq = graph.nodes[node_list[-1]]['seq']
        # seq += last_seq
        return seq


    def graph_remove_low_depth_nodes(self, must_have_st_path_after_removal=False):
        """
            Description:
                remove ANY node if depth less than cov_min_estimate
                must_have_st_path_after_removal: when true, remove node only if depth is low and there remains 1+ s-t path after removal
        """
        removed_count = 0
        flag = True
        while flag:
            flag = False
            for node in self.graph.nodes:
                if self.graph_remove_low_depth_node(node, must_have_st_path_after_removal):
                    removed_count += 1
                    flag = True
                    break
        self.logger.info("Removed {} low-depth nodes < {}.".format(removed_count, self.cov_min_estimate))
        return 0


    def graph_remove_low_depth_node(self, node, must_have_st_path_after_removal=False):
        """
            Description:
                remove the input node if depth less than cov_min_estimate
                must_have_st_path_after_removal: when true, remove node only if depth is low and there remains 1+ s-t path after removal
            Return:
                bool: whether node is removed
        """
        depth = self.graph.nodes[node]['depth'] 

        if depth >= self.cov_min_estimate:
            return False
        
        if not must_have_st_path_after_removal:
            self.logger.info('Removing node {} depth {}, length {}'.format(node, depth, str(self.graph.nodes[node]['length'])))
            self.graph.remove_node(node)
            return True
        
        # Examin if any s-t path remains after removal
        assert self.source is not None
        assert self.target is not None
        if self.source == node or self.target == node:
            return False
        if node in self.anchor_start_nodes or node in self.anchor_end_nodes:
            return False

        new_graph = self.graph.copy()
        new_graph.remove_node(node)
        if not nx.has_path(new_graph, self.source, self.target):
            self.logger.info('Keeping node {} depth {}, as removal results in no s-t path'.format(node, depth))
            return False    
        else:
            self.logger.info('Removing node {} depth {}, length {}'.format(node, depth, str(self.graph.nodes[node]['length'])))
            self.graph.remove_node(node)
            return True

    
    def assemble(self, assembler_config=None):
        """
           Description:
                A collection of main algorithms to perform assembly.
                Execute assembly with 0 mismatch in anchors. Increase number of mismatches permitted if assembly fails
                Algorithms pipeline:
                    1. identifying anchors (permitting 0 to self.max_nm_anchors mismatches)
                    2. rank anchor-pair combinations by depth
                    3. split anchors: trimming unnecessary head/tail sequences, and solve self-loop or circular path/node
                    4. find max bottleneck path + minimized path length discrepancy
        """
        flag_path_found = False
        for nm_permit in range(0, self.max_nm_anchors + 1):
            if flag_path_found:
                break
            # rank anchor combinations by depth sum
            try:
                self.get_anchor_node_lists(nm_permit)
            except nx.NodeNotFound:
                self.logger.debug("No anchor node found even if with max {} mismatches".format(nm_permit))
                self.logger.debug("Increase # mismatches allowed and try again")
                continue
            anchor_pairs = []
            for i in self.anchor_start_nodes:
                for j in self.anchor_end_nodes:
                    dp1 = self.graph.nodes[i]['depth']
                    dp2 = self.graph.nodes[j]['depth']
                    anchor_pairs.append((i, j, dp1, dp2))
            anchor_pairs.sort(key= lambda x: x[2] + x[3], reverse=True)

            # solve multiple anchor pairs problem
            for s, t, _, _ in anchor_pairs:
                self.source = s
                self.target = t
                if not nx.has_path(self.graph, self.source, self.target): 
                    continue
                
                if self.source == self.target:
                    #FIXME: need to check positions, fix trimming
                    self.contig_path = [self.source]
                    self.contig_seq = self.node_path_to_seq(self.graph, self.contig_path)
                    # return 0
                
                self.split_anchor_nodes()

                # prune graph
                self.graph_remove_low_depth_nodes(must_have_st_path_after_removal=True)
                if not nx.has_path(self.graph, self.source, self.target):   
                    continue
                
                # find optimal path
                s_t_path = self.path_max_bottleneck_weight(self.source, self.target)       
                
                self.contig_path = s_t_path
                self.contig_seq = self.node_path_to_seq(self.graph, self.contig_path)

                if s_t_path: 
                    self.logger.debug('self contig path is {}'.format(','.join(self.contig_path)))
                    self.logger.debug('self non-error-corrected & non-trimmed contig seq is {}'.format(self.contig_seq))    
                    self.logger.log("find desirable path")
                    flag_path_found = True
                    break 

                #TODO:
                # self.error_correction()
                # self.trim_contig_beyond_anchor()
                # self.logger.debug('self error-corrected and trimmed contig seq is {}'.format(self.contig_seq))
        
        return 0
    

    def path_max_bottleneck_weight(self, s, t):
        """
            Select path to maximize bottleneck weight, also least length discrepancy
            Use a min-resolution for max bottle neck, so that to tolerate some errors in bottleneck 
            
            TODO:  reimplement algorithm in https://www.biorxiv.org/content/10.1101/2021.02.26.433113v1 
                   Asofnow, the algorithm enumerates all paths, worst case is O(n!)
        """
        if not nx.has_path(self.graph, s, t):
            raise nx.NetworkXUnfeasible(
                'No path exists between s {} and t {}. Contigs does not exist in graph!'.format(s, t))
        
        min_bottle_resolution = self.cov_fav_estimate * 0.05  # covereage +- this resolution are considered equal 
        max_bottle = -1
        max_bottle_path = []
        max_bottle_length_diff = self.contig_length_fav * 1000

        for p in nx.all_simple_paths(self.graph, s, t, cutoff = self.contig_length_max):
            if self.verbose:
                self.logger.debug("checking simple s-t path p{}".format(p))
            # check inverted loop
            has_inverted_loop = False
            sorted_p = sorted(p)
            for i in range(1, len(sorted_p)):
                if sorted_p[i][:-1] == sorted_p[i-1][:-1]:
                    has_inverted_loop = True
                    break
            if has_inverted_loop:
                self.logger.debug("simple s-t path p {} has inverted loop".format(p))
                continue

            # find max bottle neck w/ min_resolution
            # for paths with the same max bottle w. min resolution, choose the path with min diff
            bottle_neck = min([self.graph.nodes[j]['depth'] for j in p])


            # not max bottle
            if bottle_neck < max_bottle - min_bottle_resolution:
                continue
            
            # examine path length 
            path_len = len(self.node_path_to_seq(self.graph, p))  
            diff = abs(path_len - self.contig_length_fav)
            if path_len < self.contig_length_min:
                self.logger.debug("simple s-t path p {} is too short, length = {}, min allowed {}".format(p, path_len, self.contig_length_min))
                continue
            if path_len > self.contig_length_max:
                self.logger.debug("simple s-t path p {} is too long, length = {}, max allowed {} ".format(p, path_len, self.contig_length_max))
                continue

            # new max bottle
            if bottle_neck > max_bottle + min_bottle_resolution:
                max_bottle = bottle_neck
                max_bottle_path = p
                max_bottle_length_diff = diff
                # print("passed bottle neck")
                continue               

            # equal w. resolution
            elif diff < max_bottle_length_diff:
                # print("same bottleneck but better length")
                max_bottle = bottle_neck
                max_bottle_path = p
                max_bottle_length_diff = diff
            else:
                # print("same bottle neck but bad length {}".format(path_len))
                pass
                    
        return max_bottle_path
    

    def write_contig(self):
        """
            only one contig is written
            contig name := name of nodes in path
        """
        output_contig_file = self.outfile
        assert self.contig_path
        if not self.contig_seq:
            self.contig_seq = self.node_path_to_seq(self.graph, self.contig_path)
        assert self.contig_seq

        with open(output_contig_file + '.fa', 'w') as f:
            trim_head = self.barcode_len if self.barcode_trim else 0            
            bc = self.contig_seq[:self.barcode_len]
            f.write('>{}_BC_{}'.format(','.join(self.contig_path), bc))
            f.write('\n')
            f.write(self.contig_seq[trim_head:])
        return 0


if __name__ == '__main__':
    # only verbose if called as main for debug purpose. Proper use will be importing and calling object
    verbose = True
    logging.basicConfig(level=logging.DEBUG)
    logging.warning('\tYou are calling AGA from main.')
    logging.warning('\tYou should NOT do so, unless for debugging.')
    logging.warning('\tYou should call via import:')
    logging.warning('\t>>> from anchorGuidedAssembler import AnchorGuidedAssembler as AGA`' + '\n' + '\n')

    parser = argparse.ArgumentParser()
    requiredarg = parser.add_argument_group('Required arguments')
    recommendarg = parser.add_argument_group('Recommended arguments')
    configarg = parser.add_argument_group('Algorithm configuration arguments')

    requiredarg.add_argument('-i', '--gfa', required=True, metavar='\b',
                                help="GFA file produced by various modes of SPAdes.")
    requiredarg.add_argument('-s1', '--anchor_start', required=True, metavar='\b')
    requiredarg.add_argument('-s2', '--anchor_end', required=True, metavar='\b')
    requiredarg.add_argument('-r1', '--read1_fq', required=True, metavar='\b')
    requiredarg.add_argument('-r2', '--read2_fq', required=True, metavar='\b')

    recommendarg.add_argument('-h1', '--ht_index_1', default="", metavar='\b',
                        help="A variable start anchor following anchor_start immediately in molecule, e.g. HT_index_1 or equivalent")
    recommendarg.add_argument('-h2', '--ht_index_2', default="", metavar='\b',
                        help="A variable end anchor following anchor_end immediately in molecule, e.g. HT_index_2 or equivalent")   
    recommendarg.add_argument('-o', '--output_contig_file', default="anchor_guide_contig", required=False, metavar='\b',
                        help='output file prefix, default: anchor_guide_contig')
    
    configarg.add_argument('--contig_barcode_len', type=int, metavar='\b',
                        help='length of contig barcode')
    
    configarg.add_argument('--trim_barcode', type=bool, metavar='\b',
                        help='if trimming barcode from contig ')

    configarg.add_argument('--max_nm_anchors', type=int, default=2, metavar='\b',
                        help='maximum number of anchors permitted')

    configarg.add_argument('--use_ambiguous_anchor', type=bool, metavar='\b',
                        help='if anchor node is ambiguous, use the max-depth one')
    configarg.add_argument('--allow_repeats_same_orientation', type=bool, metavar='\b')
    configarg.add_argument('--allow_repeats_revcomp', type=bool, metavar='\b')
    configarg.add_argument('--always_keep_ndoe_longer_than_length_good', type=bool, metavar='\b')
    
    args = parser.parse_args()


    print(args.read1_fq)
    left= args.read1_fq
    AGA = AnchorGuidedAssembler(
                        args.gfa,
                        args.anchor_start,
                        args.anchor_end,
                        left_fq=args.read1_fq,
                        right_fq=args.read2_fq,
                        output_contig_file = args.output_contig_file,
                        ht_index_1=args.ht_index_1,
                        ht_index_2=args.ht_index_2,
                        barcode_trim = True,
                        max_nm_anchors = args.max_nm_anchors,
                        barcode_len = args.contig_barcode_len,   
                        verbose=verbose)