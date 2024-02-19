#! /usr/bin/env python3

import os
import shutil
import re
import argparse
import random
import pysam


R1_ORIENTATION_FORWARD = 'forward'
R1_ORIENTATION_REVERSE = 'reverse'
R1_ORIENTATION_UNKNOWN = 'unknown'

def lcs(R, T, i, j, mismatch=0, mismatch_threshold=1):
    if mismatch > mismatch_threshold:
        return 0
    elif i == len(R) or j == len(T):
        return 0
    elif re.match(R[i], T[j]):
        return 1 + lcs(R, T, i + 1, j + 1, mismatch, mismatch_threshold)
    else:
        mismatch += 1
        return max(lcs(R, T, i, j + 1, mismatch, mismatch_threshold),
                   lcs(R, T, i + 1, j, mismatch, mismatch_threshold),
                   lcs(R, T, i + 1, j + 1, mismatch, mismatch_threshold))

class IupacPattern:
    IUPAC = {'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT', 'V': 'ACG'}

    def __init__(self, pattern_list=None, pattern=None, rc=False):
        if pattern_list is None:
            if pattern is None:
                raise ValueError('One of pattern or pattern_list must be supplied')
            pattern_list = [pattern]
        if rc:
            pattern_list = [revcomp(pat) for pat in pattern_list]
        self.pattern_list = pattern_list
        self.iupac_patterns = [self.iupac_pattern(pat) for pat in pattern_list]

    @staticmethod
    def iupac(code):
        return IupacPattern.IUPAC[code]

    @staticmethod
    def iupac_pattern(pattern):
        return [r'{}'.format(c) if c in 'ATGCN' else '[{}]'.format(IupacPattern.iupac(c)) for c in pattern]

    def find_pattern(self, seq, threshold=0, five_prime=True, search_start=None, search_end=None):
        for which, pattern in enumerate(self.iupac_patterns):
            if len(seq) > len(pattern):
                max_score = 0
                start = -1
                stop = -1
                if search_start is None:
                    if five_prime:
                        # start at 0, go 1/2 way
                        search_start = 0
                        search_end = int(len(seq) / 2) - len(pattern) + 1
                    else:
                        # start 1/2 way thru, go to end
                        search_start = int(len(seq) / 2)
                        search_end = len(seq) - len(pattern) + 1
                elif search_end is None:
                    raise ValueError('if search_start is supplied, search_end must also be supplied')

                for i in range(search_start, search_end):
                    subseq = seq[i:i + len(pattern)]
                    score = lcs(pattern, subseq, 0, 0, mismatch_threshold=threshold)
                    # keep the last max if coming from 5',
                    # first max if coming from 1/2 thru to the end
                    if (five_prime and score >= max_score) or \
                            (not five_prime and score > max_score):
                        max_score = score
                        start = i
                        stop = i + len(pattern)

                min_score = len(pattern) - threshold
                if max_score >= min_score:
                    return start, stop, which, True

        return -1, -1, -1, False

    def trim_anchor_start(self, seq, mismatch_threshold=1, anchor_remove=False):
        start, stop, matched, found_anchor = self.find_pattern(seq, threshold=mismatch_threshold, five_prime=True)
        if found_anchor:
            if anchor_remove:
                return seq[stop:], seq[0:stop], matched, True
            else:
                return seq[start:], seq[0:start], matched, True
        return seq, '', -1, False

    def trim_anchor_end(self, seq, mismatch_threshold=1, anchor_remove=False):
        start, stop, matched, found_anchor = self.find_pattern(seq, threshold=mismatch_threshold, five_prime=False)
        if found_anchor:
            if anchor_remove:
                return seq[:start], seq[start:], matched, True
            else:
                return seq[:stop], seq[stop:], matched, True
        return seq, '', -1, False

    def find_anchor_prefix(self, seq, start_from=0, s_range=60, mismatch_threshold=0, length=12):
        start, stop, matched, found_anchor = self.find_pattern(
            seq, threshold=mismatch_threshold, five_prime=True, search_start=start_from, search_end=start_from + s_range
        )
        if found_anchor:
            return seq[start - length:start]
        else:
            return ''

    def find_anchor(self, seq, start_from=0, s_range=60, mismatch_threshold=0, length=12):
        start, stop, matched, found_anchor = self.find_pattern(
            seq, threshold=mismatch_threshold, five_prime=True, search_start=start_from, search_end=start_from + s_range
        )
        if found_anchor:
            return seq[stop:stop + length]
        else:
            return ''

dna_mapping = str.maketrans('ACGTNRYSWKMBV', 'TGCANYRSWMKAT')

def get_strand_from_short_read(contig, l_file, tmp_dir, threads=1):
    fafile = f"{tmp_dir}/fafile.fa"
    with open(fafile, "w") as output_handle:
        output_handle.write(">ref\n")
        output_handle.write(contig + "\n")

    command = f"bwa-mem2 index {fafile}"
    retcode = os.system(command)

    if retcode != 0:
        return R1_ORIENTATION_UNKNOWN

    command = f"bwa-mem2 mem -t {threads} {fafile} {l_file} > {tmp_dir}/output.sam"
    retcode = os.system(command)
    if retcode != 0:
        return R1_ORIENTATION_UNKNOWN

    aligned_count = 0
    forward_count = 0
    for segment in pysam.AlignmentFile(f"{tmp_dir}/output.sam"):
        if not segment.is_unmapped:
            if not (segment.is_supplementary or segment.is_secondary):
                aligned_count += 1
                if not segment.is_reverse:
                    forward_count += 1
    
    if aligned_count == 0:
        return R1_ORIENTATION_UNKNOWN
    
    if forward_count / float(aligned_count) > 0.5:
        return R1_ORIENTATION_FORWARD
    else:
        return R1_ORIENTATION_REVERSE


def revcomp(seq, qual=None):
    """

    :param seq:
    :param qual:
    :return:
    """
    if isinstance(seq, str):
        seq = seq.translate(dna_mapping)[::-1]
    elif isinstance(seq, list):
        seq = [x.translate(dna_mapping) for x in seq][::-1]
    if qual is not None:
        qual = qual[::-1]
        return seq, qual
    return seq

def trim_contig(contig, term_len, r1_frag, pcr_primer_list):
    term_len = min(term_len, int(len(contig) / 2))

    three_term_seq = contig[-term_len:]
    five_term_seq = contig[0:term_len]
    middle_seq = contig[term_len:len(contig) - term_len]

    len_original_3prime = len(three_term_seq)
    len_original_5prime = len(five_term_seq)

    # Check for r1_frag in 3' tail and if present, trim it and everything after
    index = three_term_seq.find(r1_frag)
    if index != -1:
        three_term_seq = three_term_seq[:index]
    else:
        # If not present, check for r1_frag_rc in 5' and trim if found
        r1_frag_rc = revcomp(r1_frag)
        index = five_term_seq.find(r1_frag_rc)
        if index != -1:
            five_term_seq = five_term_seq[index + len(r1_frag_rc):]
        else:
            # check for reverse orientation
            index = five_term_seq.find(r1_frag)
            if index != -1:
                five_term_seq = five_term_seq[index + len(r1_frag):]
            else:
                index = three_term_seq.find(r1_frag_rc)
                if index != -1:
                    three_term_seq = three_term_seq[:index]

    # Check for last 9 bases of the PCR primer from config in 5' terminal
    pcr_primer = ''
    pcr_primer_rc = ''
    if len(pcr_primer_list) == 1:
        pcr_primer = pcr_primer_list[0]
        pcr_primer_rc = revcomp(pcr_primer_list[0])
    else:
        pcr_primer = pcr_primer_list[1]
        pcr_primer_rc = revcomp(pcr_primer_list[0])

    pcr_frag = pcr_primer[-7:]
    index = five_term_seq.find(pcr_frag)
    if index != -1:
        five_term_seq = five_term_seq[index + len(pcr_frag):]
    else:
        # If not present, check for the pcr_frag_rc in the 3' terminal and trim it
        pcr_frag_rc = revcomp(pcr_frag)
        index = three_term_seq.find(pcr_frag_rc)
        if index != -1:
            three_term_seq = three_term_seq[:index]
        else:
            # check for reverse orientation
            index = three_term_seq.find(pcr_frag)
            if index != -1:
                # if present, extend the match for as long as we can until hitting a mismatch.
                trim_len = len(pcr_frag)
                while three_term_seq[index - 1:].startswith(pcr_primer[-trim_len - 1:]) and trim_len < len(pcr_primer):
                    index -= 1
                    trim_len += 1
                # print('found fragment {} at index {}'.format(pcr_primer[-trim_len:], index))
                three_term_seq = three_term_seq[:index]
            else:
                index = five_term_seq.find(pcr_frag_rc)
                if index != -1:
                    trim_len = len(pcr_frag_rc)
                    while five_term_seq[:index + trim_len + 1].endswith(pcr_primer_rc[:trim_len + 1]) and trim_len < len(pcr_primer_rc):
                        trim_len += 1
                    # print('found fragment {} at index {}'.format(pcr_primer_rc[:trim_len], index))
                    five_term_seq = five_term_seq[index + trim_len:]

    contig_trimmed = five_term_seq + middle_seq + three_term_seq
    len_5prime_trimmed = len_original_5prime - len(five_term_seq)
    len_3prime_trimmed = len_original_3prime - len(three_term_seq)

    return contig_trimmed, len_5prime_trimmed, len_3prime_trimmed

def trim_anchor_seqs(seq, anchor_start, anchor_end, anchor_remove=True):  
    trimmed_5prime_seq = ''
    trimmed_3prime_seq = ''
    is5pDetected = False
    is3pDetected = False
    anchor_start_idx = 0
    anchor_end_idx = 0

    if anchor_start:
        matcher = IupacPattern(pattern_list=[anchor_start,], rc=False)
        seq, trimmed_5prime_seq, anchor_start_idx, is5pDetected = matcher.trim_anchor_start(
            seq, anchor_remove=anchor_remove
        )

    if anchor_end:
        matcher = IupacPattern(pattern_list=[anchor_end,], rc=False)
        seq, trimmed_3prime_seq, anchor_end_idx, is3pDetected = matcher.trim_anchor_end(
            seq, anchor_remove=anchor_remove
        )

    return seq, trimmed_5prime_seq, trimmed_3prime_seq, is5pDetected, is3pDetected, anchor_start_idx, anchor_end_idx

def trim_all_anchors(seq, anchor_start, anchor_end, ht_anchor_start):
    seq, _, _, is5pDetected, is3pDetected, _, _ = trim_anchor_seqs(seq, anchor_start, anchor_end)
    if is5pDetected and ht_anchor_start is not None and seq.startswith(ht_anchor_start):
        seq = seq[len(ht_anchor_start):]
    return seq, is5pDetected, is3pDetected

def read_fasta(input_filename):
    result = []
    current_entry = ""
    for line in open(input_filename):
        if line.startswith(">"):
            if len(current_entry) > 0:
                result.append(current_entry)
            current_entry = ""
        else:
            current_entry += line.rstrip()
    if len(current_entry) > 0:
        result.append(current_entry)
    return result

def status(trim5, trim3):
    if trim5 and trim3:
        return "Full-length"
    elif trim5 and not trim3:
        return "Start_only"
    elif not trim5 and trim3:
        return "End_only"
    else:
        return "Undetected"

def rank_status(status):
    if status == "Full-length":
        return 2
    if status == "Undetected":
        return 0
    return 1

def update_final_contig(final_contig, final_status, current_contig, current_status):
    final_rank = rank_status(final_status)
    current_rank = rank_status(current_status)

    if current_rank > final_rank:
        return current_contig, current_status
    elif final_rank > current_rank:
        return final_contig, final_status
    else:
        if len(final_contig) > len(current_contig):
            return final_contig, final_status
        return current_contig, current_status

def read_fastq(input_file):
    lines = []
    for line in open(input_file):
        lines.append(line.rstrip())

    result = []
    for idx in range(int(len(lines)/4)):
        result.append((lines[4*idx], lines[4*idx + 1], lines[4*idx + 2], lines[4*idx + 3]))

    return result

def write_fastq(entries, output_file):
    with open(output_file, "w") as output_handle:
        for entry in entries:
            for element in entry:
                output_handle.write(element + "\n")

def sample_fastq(input_file_forward, input_file_reverse, sampling_target, output_file_forward, output_file_reverse):
    forward_entries = read_fastq(input_file_forward)
    reverse_entries = read_fastq(input_file_reverse)
    if sampling_target >= len(forward_entries):
        return (input_file_forward, input_file_reverse)
    
    final_entries = random.sample(list([entry for entry in zip(forward_entries, reverse_entries)]), sampling_target)
    os.makedirs(os.path.dirname(output_file_forward), exist_ok=True)
    os.makedirs(os.path.dirname(output_file_reverse), exist_ok=True)
    write_fastq([entry[0] for entry in final_entries], output_file_forward)
    write_fastq([entry[1] for entry in final_entries], output_file_reverse)
    return output_file_forward, output_file_reverse


def do_assembly_iteration(forward_fq, reverse_fq, threads, sampling_target, pcr_primer, r1_orientation):
    if sampling_target is not None:
        forward_fq, reverse_fq = sample_fastq(forward_fq, reverse_fq, sampling_target, f"spades_output/forward.fq", f"spades_output/reverse.fq")

    #
    # run assembly with spades
    #
    command = f"spades.py -k 21,33,55,77,99,127 -t {threads} --careful --sc -o spades_output --phred-offset 33 --disable-gzip-output -1 {forward_fq} -2 {reverse_fq}"
    retcode = os.system(command) 

    #
    # read and trim entries
    #
    trimmed_contigs = []
    if retcode == 0 and os.path.isfile("spades_output/contigs.fasta"):
        for current_contig in read_fasta("spades_output/contigs.fasta"):
            current_contig, _, _ = trim_contig(current_contig, 40, "CCTACAC", pcr_primer)
            trimmed_contigs.append(current_contig)
    
    if r1_orientation != R1_ORIENTATION_UNKNOWN and len(trimmed_contigs) > 0:
        contig_orientation = get_strand_from_short_read(trimmed_contigs[0], forward_fq, "spades_output", threads)
        if contig_orientation != R1_ORIENTATION_UNKNOWN:
            if contig_orientation != r1_orientation:
                trimmed_contigs = list([revcomp(trimmed_contig) for trimmed_contig in trimmed_contigs])

    shutil.rmtree("spades_output", ignore_errors=True)
    return trimmed_contigs

def create_solo_barcodes_dict(solo_barcodes, solo_contig_barcodes):
    result = {}
    if solo_barcodes is not None and len(solo_barcodes) > 0:
        for solo_barcode, solo_contig_barcode in zip(solo_barcodes.rstrip().split(","), solo_contig_barcodes.rstrip().split(",")):
            result[solo_barcode] = solo_contig_barcode
    return result


def spades_assembly(input_dir, output_prefix, pcr_primer, anchor_start, anchor_end, threads, assembly_iterations, sampling_target, solo_barcodes, solo_contig_barcodes, r1_orientation):
    solo_barcodes_dict = create_solo_barcodes_dict(solo_barcodes, solo_contig_barcodes)    
    with open(f"{output_prefix}_output.fa", "w") as output_handle, open(f"{output_prefix}_output.csv", "w") as csv_handle:
        for candidate_file in os.listdir(input_dir):
            if candidate_file.endswith("_R1.fastq"):
                forward_fq = candidate_file
                reverse_fq = candidate_file.replace("_R1.fastq", "_R2.fastq")
                umi = candidate_file.split("_")[0]

                if len(solo_barcodes) > 0:
                    if umi not in solo_barcodes_dict:
                        continue
                    else:
                        contig_barcode = solo_barcodes_dict[umi]
                else:
                    contig_barcode = None
                
                input_forward_fq = f"{input_dir}/{forward_fq}"
                input_reverse_fq = f"{input_dir}/{reverse_fq}"

                read_count = len(read_fastq(input_forward_fq))
                #
                # do a number of assembly iterations and sort the contigs
                # by number of observations across iterations
                #
                contigs = []
                contig_counts = {}
                for _ in range(assembly_iterations if sampling_target is not None and sampling_target < read_count else 1):
                    for contig in do_assembly_iteration(input_forward_fq, input_reverse_fq, threads, sampling_target, pcr_primer, r1_orientation):
                        if contig in contig_counts:
                            contig_counts[contig] += 1
                        else:
                            contigs.append(contig)
                            contig_counts[contig] = 1
                contigs.sort(key = lambda x: -contig_counts[x])

                #
                # Find the most frequent Full-length contig
                #
                (final_contig, final_status) = (None, None)
                (first_contig, first_status) = (None, None)
                for current_contig in contigs:
                    current_contig, trim5, trim3 = trim_all_anchors(current_contig, anchor_start if len(anchor_start) > 0 else None, anchor_end if len(anchor_end) > 0 else None, contig_barcode)
                    current_status = status(trim5, trim3)
                    if first_contig is None:
                        (first_contig, first_status) = (current_contig, current_status)
                    if current_status == "Full-length":
                        (final_contig, final_status) = (current_contig, current_status)
                        break
                    if current_status != "Undetected" and final_contig is None:
                        (final_contig, final_status) = (current_contig, current_status)

                if final_contig is None:
                    (final_contig, final_status) = (first_contig, first_status)
                
                if final_contig is not None:
                    output_handle.write(f">{umi}_{len(final_contig)}\n")
                    output_handle.write(final_contig + "\n")
                    if sampling_target is not None and sampling_target < read_count:
                        effective_read_count = sampling_target
                    else:
                        effective_read_count = read_count
                    csv_handle.write(f"{umi},{len(final_contig)},{effective_read_count},{final_status}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Assemble LoopSeq UMI fragments with Spades")
    parser.add_argument("--pcr-primer", dest="pcr_primer", required=True, help="Nucleotide sequence of PCR primer")
    parser.add_argument("--anchor-start", dest="anchor_start", required=True, help="Expected start sequence of long read")
    parser.add_argument("--anchor-end", dest="anchor_end", required=True, help="Expected end sequence of long read")
    parser.add_argument("--threads", type = int, dest="threads", default = 1, help="Number of threads")
    parser.add_argument("--assembly-iterations", type=int, dest="assembly_iterations", default = 1, help="Number of iterations to repeat assembly (for Solo)")
    parser.add_argument("--sampling-target", type=int, dest="sampling_target", default=None, help="Down-sampling target for each assembly iteration (for Solo)")
    parser.add_argument("--solo-contig-barcodes", dest="solo_contig_barcodes", default="", help="Comma-delimited list of contig barcodes (for Solo)")
    parser.add_argument("--solo-barcodes", dest="solo_barcodes", default="", help="Comma-delimited list of Solo barcodes")
    parser.add_argument("--r1-orientation", dest="r1_orientation", default = R1_ORIENTATION_UNKNOWN, help="Orientation of R1 short reads (if known), relative to long read (FORWARD, REVERSE, or UNKNOWN)")
    parser.add_argument("input_dir")
    parser.add_argument("output_prefix")

    args = parser.parse_args()

    spades_assembly(args.input_dir, args.output_prefix, args.pcr_primer, args.anchor_start, args.anchor_end, args.threads, args.assembly_iterations, args.sampling_target, args.solo_barcodes, args.solo_contig_barcodes, args.r1_orientation)
