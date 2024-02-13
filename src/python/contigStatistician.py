"""
Author: Carl Zang
Date Started: Jul 2023
Last updated: Aug 1 2023 or check git history 

Version:    v1.0.0

Compute k-mer coverage. Estimate contig coverage and length from LoopSeq reads.

Usage: see test() for an example. Import module/class. Don't run from main.


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


from loop.common.utils import RunShellCommand
from collections import Counter
import sys
import os


class FqStatistician:
    def __init__(self, fqfile, nameprefix="", k=33):
        # config param
        self.fqfile = fqfile
        self.nameprefix = nameprefix + "FqStat_Estimation"
        self.kmc_file = ""
        self.k = k
        # Stats
        self.kmers_depth = []
        self.depth_counter = None
        self.num_seqs = -1
        self.num_bases = -1
        self.read_avg_len = -1
        self.read_min_len = -1
        self.read_max_len = -1
        self.effective_len = -1
        self.get_kmer_depth()
        self.get_depth_counter()
        self.get_fq_stats()
        # Estimates
        self.cov_lower_bound = -1
        self.cov_upper_bound = -1
        self.cov_estimation  = -1
        self.cov_50percentile_cut_tail = -1
        self.estimate_coverage()
    

    def get_kmer_depth(self):
        """
            Call kmc to compute kmer depth
        """
        fqfile = self.fqfile
        if fqfile.endswith('.fq') or fqfile.endswith('.fastq'):
            fqflag = '-fq {}'.format(fqfile)
        elif fqfile.endswith('.fa') or fqfile.endswith('.fasta'):
            fqflag = '-fa {}'.format(fqfile)
        else:
            raise IOError("File {} is not fasta nor fastq file!".format(self.fqfile))
        
        # build KMC database
        kmc_executable = 'kmc'
        k = self.k
        kflag = '-k' + str(k)
        exclude_less_than = '-ci0'
        max_value_of_counter = '-cs100000000' #1e8
        donot_use_canonical = '-b'
        fqbase = os.path.basename(fqfile)
        output_kmc_db = '{}.{}.k{}'.format(self.nameprefix, fqbase, k)
        work_dir = "."
        cmd_kmcbuild = ' '.join([kmc_executable, kflag, exclude_less_than, max_value_of_counter, donot_use_canonical, fqflag, output_kmc_db, work_dir])
        RunShellCommand(cmd_kmcbuild, 'Build KMC database for {} {} {} {}'.format(fqfile, kflag, exclude_less_than, max_value_of_counter))

        # dump and sort KMC database
        kmc_dump_executable = 'kmc_dump'
        output_kmc_dumped = '{}.{}.k{}.kmc.tsv'.format(self.nameprefix, fqbase, k)
        cmd_kmcdump = ' '.join([kmc_dump_executable, exclude_less_than, output_kmc_db, output_kmc_dumped])
        RunShellCommand(cmd_kmcdump, 'Dump KMC database for {} {} {} {}'.format(fqfile, k, exclude_less_than, max_value_of_counter))

        output_kmc_sorted = '{}.{}.k{}.sorted.tsv'.format(self.nameprefix, fqbase, k)
        cmd_kmc_sort = 'sort -r -k2 -n {} > {}'.format(output_kmc_dumped, output_kmc_sorted)
        RunShellCommand(cmd_kmc_sort, 'Sorting KMC output')

        # Read
        with open(output_kmc_sorted, 'r') as f:
            for l in f.readlines():
                kmer, count = l.split()
                count = int(count)
                self.kmers_depth.append((kmer, count))
        
        assert len(self.kmers_depth) > 0
        for i in range(1, len(self.kmers_depth)):
            assert self.kmers_depth[i][1] <= self.kmers_depth[i-1][1]

        return 0
    
    def get_depth_counter(self):
        self.depth_counter = Counter([x[1] for x in self.kmers_depth])

    def get_fq_stats(self):
        seqkit_stat = 'seqkit stat {}'.format(self.fqfile)
        statline = RunShellCommand(seqkit_stat, "Get {} stats".format(self.fqfile), return_stdout=True).split('\n')[1].split()
        #file     format  type  num_seqs  sum_len  min_len  avg_len  max_len
        self.num_seqs = int(statline[3].replace(',', ''))
        self.num_bases = int(statline[4].replace(',', ''))
        self.read_avg_len = float(statline[6].replace(',', ''))
        self.read_min_len = float(statline[5].replace(',', ''))
        self.read_max_len = float(statline[7].replace(',', ''))
        self.effective_len = self.read_avg_len - self.k + 1
        return 0

    def n_percentage_cover_kmer(self, percentage = 0.5):
        """
            return the (kmer,cov) for those kmers with coverage higher than (incld.) this kmer accounts for more than %percentage of all kmers
        """
        total_cov = sum([x[1] for x in self.kmers_depth])
        current_cov = 0
        assert percentage >= 0 and percentage <= 1
        percentage_cov = total_cov * percentage
        for kmer, cov in self.kmers_depth:
            current_cov += cov
            if current_cov >= percentage_cov:
                return kmer,cov
        raise ValueError("Something ungodly wrong of the math happened!")
    
    def estimate_coverage(self):
        """
            Description: estimation is based on kmer cov
                incld. n50 kmer, 50 percentile kmer
            Note: kmer frequency has very long tail. Using n25 is much better than 25 percentile
        """
        n25_kmer = self.n_percentage_cover_kmer(0.25)
        n50_kmer = self.n_percentage_cover_kmer(0.5)
        n75_kmer = self.n_percentage_cover_kmer(0.75)
        percentile25_kmer = self.kmers_depth[len(self.kmers_depth)//4]
        percentile75_kmer = self.kmers_depth[len(self.kmers_depth)*3//4]

        # trade off percentile and n_percentage
        self.cov_lower_bound = min(n75_kmer[1], percentile75_kmer[1]) 
        self.cov_upper_bound = max(n25_kmer[1], percentile25_kmer[1])
        self.cov_estimation = n50_kmer[1]

        n95_kmer = self.n_percentage_cover_kmer(0.80)
        num_n95_kmer = [x for x in self.kmers_depth if x[1] > n95_kmer[1]]
        self.cov_50percentile_cut_tail = num_n95_kmer[len(num_n95_kmer)//2][1]

        if self.cov_estimation < 5:
            print("Estimated Coverage {} is too low".format(self.cov_estimation))
        return 0


class ContigStatistician:
    def __init__(self, left, right, nameprefix="", k=33):
        """
            Estimate contig statistics from reads (spacer & HT_index_2 removed from R1)
            Parameters:
                k: kmer size
        """
        # stats
        self.left_stats = FqStatistician(left, nameprefix, k)
        self.right_stats = FqStatistician(right, nameprefix, k)
        self.k = k

        # Estimates
        self.nameprefix = nameprefix + "ContigStat_Estimation"
        self.cov_upper_bound = 0
        self.cov_lower_bound = 0
        self.cov_estimation = 0
        self.cov_50percentile_cut_tail = 0
        self.num_bases = 0
        self.length_upper_bound = 0
        self.length_lower_bound = 0
        self.length_estimation = 0
        self.length_estimation_by_cut_tail = 0
        self.length_good = 0
        self.estimate_contig_coverage()
        self.estimate_contig_length_range()
        self.leverage_estimetes()
        self.write_stats()

    def contig_range(self):
        return self.length_good, self.length_lower_bound, self.length_upper_bound
    
    def contig_length(self):
        return self.length_good 
    
    def estimate_contig_coverage(self):
        """
            Description: estimation is based on kmer cov
                incld. n50 kmer, 50 percentile kmer
        """
        for fqstats in [self.left_stats, self.right_stats]:
            coeff = fqstats.read_avg_len / fqstats.effective_len
            self.cov_upper_bound += fqstats.cov_upper_bound * coeff
            self.cov_lower_bound += fqstats.cov_lower_bound * coeff
            self.cov_estimation  += fqstats.cov_estimation  * coeff 
            self.cov_50percentile_cut_tail += fqstats.cov_50percentile_cut_tail * coeff

        assert self.cov_lower_bound >= 0
        assert self.cov_lower_bound <= self.cov_estimation
        assert self.cov_upper_bound >= self.cov_estimation

        return 0

    def estimate_contig_length_range(self):
        """
            estimate: lower/upper bound and good estimation of contig length range 
        """
        self.num_bases = self.left_stats.num_bases + self.right_stats.num_bases 

        self.length_lower_bound = self.num_bases / self.cov_upper_bound
        self.length_upper_bound = self.num_bases / self.cov_lower_bound
        self.length_estimation = self.num_bases / self.cov_estimation
        self.length_estimation_by_cut_tail = self.num_bases / self.cov_50percentile_cut_tail

        if self.length_estimation < 500:
            print("Estimated length {} from n50 kmer is too short".format(self.length_estimation))
        if self.length_estimation > 10000:
            print("Estimated length {} from n50 kmer is too long".format(self.length_estimation))
        if self.length_estimation_by_cut_tail < 500:
            print("Estimated length {} from cutting tail is too short".format(self.length_estimation_by_cut_tail))
        if self.length_estimation_by_cut_tail > 10000:
            print("Estimated length {} from cutting tail is too long".format(self.length_estimation_by_cut_tail))
        return 0    
    
    def leverage_estimetes(self):
        """
            Some aggressive heuristics to set upper/lower bound more tight
        """
        # length_estimation_by_cut_tail is heavily affected by low frequency k-mers. Error rate can go as high as 35%. Do not use.

        # is_length_estimation_bad             = self.length_estimation < 500             or self.length_estimation > 10000
        # is_length_estimation_by_cut_tail_bad = self.length_estimation_by_cut_tail < 500 or self.length_estimation_by_cut_tail > 10000
        # if is_length_estimation_bad       and not is_length_estimation_by_cut_tail_bad:
        #     self.length_good = self.length_estimation_by_cut_tail
        # elif not is_length_estimation_bad and     is_length_estimation_by_cut_tail_bad:
        #     self.length_good = self.length_estimation
        # else:
        #     self.length_good = (self.length_estimation_by_cut_tail + self.length_estimation) / 2
        self.length_good = self.length_estimation
        self.length_lower_bound = min(max(self.length_lower_bound, 0.5 * self.length_good), self.length_estimation * 0.75)#, self.length_estimation_by_cut_tail * 0.75)
        self.length_upper_bound = max(min(self.length_upper_bound, 2   * self.length_good), self.length_estimation * 1.25)#, self.length_estimation_by_cut_tail * 1.25)
        return 0

    def write_stats(self):
        n = self.nameprefix + ".stats.tsv"
        with open(n, 'w') as f:
            l1 = '{}\t{}\t{}\t{}\t{}\t{}'.format('name','estimate', 'median_no_tail', 'lower_bound', 'upper_bound', 'recommend')
            l2 = '{}\t{}\t{}\t{}\t{}\t{}'.format('left.cov',    self.left_stats.cov_estimation,     self.left_stats.cov_50percentile_cut_tail,  self.left_stats.cov_lower_bound,    self.left_stats.cov_upper_bound,        'NA')
            l3 = '{}\t{}\t{}\t{}\t{}\t{}'.format('right.cov',   self.right_stats.cov_estimation,    self.right_stats.cov_50percentile_cut_tail, self.right_stats.cov_lower_bound,   self.right_stats.cov_upper_bound,       'NA')
            l4 = '{}\t{}\t{}\t{}\t{}\t{}'.format('contig.cov',  self.cov_estimation,                self.cov_50percentile_cut_tail,             self.cov_lower_bound,               self.cov_upper_bound,                   'NA')
            l5 = '{}\t{}\t{}\t{}\t{}\t{}'.format('contig.len',  self.length_estimation,             self.length_estimation_by_cut_tail,         self.length_upper_bound,            self.length_lower_bound,                self.length_good)
            f.write('\n'.join([l1, l2, l3, l4, l5]))

def test():
    # Usage:
    #   cs = ContigStatistician(left_path, right_path, 'out_name_prefix', k = 33)
    #   print("estimated contig length is approx {}".format(cs.length_good))
    #
    # for debug: python contigStatistician.py test [left.fq] [right.fq]
    
    if len(sys.argv) == 4:
        left = sys.argv[2] 
        right = sys.argv[3]
    elif len(sys.argv) == 2:
        left = 'scratch_diff_read_coverage/left.500.fq'
        right = 'scratch_diff_read_coverage/right.500.fq'
    else:
        print("wrong number of arguments!", file = sys.stderr)

    cs = ContigStatistician(left, right, 'test', k = 33)
    print('{}\t{}\t{}\t{}\t{}\t{}'.format('name','est', '50pct_no_tail', 'lower', 'upper', 'good'))
    print('{}\t{}\t{}\t{}\t{}\t{}'.format('R1.cov', cs.left_stats.cov_estimation, cs.left_stats.cov_50percentile_cut_tail, cs.left_stats.cov_lower_bound, cs.left_stats.cov_upper_bound, 'NA'))
    print('{}\t{}\t{}\t{}\t{}\t{}'.format('R2.cov', cs.right_stats.cov_estimation, cs.right_stats.cov_50percentile_cut_tail, cs.right_stats.cov_lower_bound, cs.right_stats.cov_upper_bound, 'NA'))
    print('{}\t{}\t{}\t{}\t{}\t{}'.format('cs.cov', cs.cov_estimation, cs.cov_50percentile_cut_tail, cs.cov_lower_bound, cs.cov_upper_bound, 'NA'))
    print('{}\t{}\t{}\t{}\t{}\t{}'.format('cs.len',cs.length_estimation, cs.length_estimation_by_cut_tail, cs.length_upper_bound, cs.length_lower_bound, cs.length_good))
    print("ContigStatistician test finished!")

if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] in ['test', '--test', '-test']:
        print("Debug tests for contigStatistician:")
        test()
    else:
        raise NotImplementedError("You should not run ContigStatistician from main. Import and call the object ")