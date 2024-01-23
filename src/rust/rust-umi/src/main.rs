use anyhow::{Ok, Result};
use clap::Parser;
use itertools::Itertools;
use std::cmp;
use std::collections::HashMap;
use std::fs::{create_dir, File};
use std::io::{self, BufRead, BufWriter, Write};
use std::iter::zip;
use std::path::Path;

#[cfg(test)]
mod tests {
    use crate::trim_forward_read;
    use crate::trim_reverse_read;

    #[test]
    fn test_trim_forward_read() {
        //
        // Test simple matching
        //
        assert_eq!(
            &"ACGTAAAA"[trim_forward_read("ACGTAAAA", "GT", 2)..],
            "AAAA"
        );

        //
        // Test partial match at the end of the string
        //
        assert_eq!(
            &"ACGTATT"[trim_forward_read("ACGTATT", "GTATTATT", 2)..],
            ""
        );
    }

    ///trim_reverse_read(reverse_sequence : &str, trim_sequence : &str) -> usize
    #[test]
    fn test_trim_reverse_read() {
        //
        // Test simple matching
        //
        assert_eq!(
            &"AAAAAGTGTGAAAA"[..trim_reverse_read("AAAAAGTGTGAAAA", "GTGTG")],
            "AAAAA"
        );

        //
        // Test too short for match
        //
        assert_eq!(
            &"AAAAAGTGTGAAAA"[..trim_reverse_read("AAAAAGTGTGAAAA", "GTGT")],
            "AAAAAGTGTGAAAA"
        );
        assert_eq!(&"ACTG"[..trim_reverse_read("ACTG", "ACTGACTG")], "ACTG");

        //
        // Test seed match failed
        //
        assert_eq!(
            &"AAAAAGTGTGTTTAAAA"[..trim_reverse_read("AAAAAGTGTGTTTAAAA", "GTGTGTCC")],
            "AAAAAGTGTGTTTAAAA"
        );

        //
        // Test partial match
        //
        assert_eq!(
            &"AAAAAGTGTGTTTGTAAAA"[..trim_reverse_read("AAAAAGTGTGTTTGTAAAA", "GTGTGTTCAT")],
            "AAAAA"
        );
    }
}

#[derive(Parser, Debug)]
#[clap(author="Ryan Kelley", version, about="Trim and group by UMI", long_about = None)]
struct Args {
    /// Uncompresed FASTQ file with paired forward reads
    forward_paired: String,

    /// Uncompresed FASTQ file with paired reverse reads
    reverse_paired: String,

    /// Uncompresed FASTQ file with unpaired forward reads (currently ignored)
    forward_unpaired: String,

    // assigned sample index
    sample_index: String,

    /// Nucleotide sequence of PCR primer
    #[clap(short, long, default_value = "")]
    pcr_primer: String,

    /// Starting position to search for barcode
    #[clap(short, long, default_value_t = 16)]
    umi_length: usize,

    /// Minimum length of trimmed sequence to output
    #[clap(long, default_value_t = 30)]
    min_length: usize,

    /// Minimum number of sequences to output per UMI
    #[clap(long, default_value_t = 100)]
    min_count: usize,

    /// Enable processing of Solo data
    #[clap(long, takes_value = false)]
    solo: bool,

    /// Nucleotide sequence of barcode spacer
    #[clap(long)]
    barcode_spacer: Option<String>,
}

fn switch_base(c: char) -> char {
    match c {
        'A' => 'T',
        'C' => 'G',
        'T' => 'A',
        'G' => 'C',
        'U' => 'A',
        _ => 'N',
    }
}

fn find_solo_plate_index_position<'a, 'b>(seq: &'a str, bc_len: usize, spacer: &'b str) -> usize {
    for idx in 0..4 {
        let spacer_start = bc_len + idx + 1;
        let spacer_end = spacer_start + spacer.len();
        if spacer_end < seq.len() && &seq[spacer_start..spacer_end] == spacer {
            return spacer_end;
        }
    }
    return 0;
}

fn revcomp(dna: &str) -> String {
    // result vector
    let mut rdna: String = String::with_capacity(dna.len());

    // iterate through the input &str
    for c in dna.chars().rev() {
        rdna.push(switch_base(c));
    }
    return rdna;
}

fn trim_forward_read(
    forward_sequence: &str,
    rc_primer_sequence: &str,
    search_start: usize,
) -> usize {
    let max_match_count = cmp::min(
        rc_primer_sequence.len(),
        forward_sequence.len() - search_start,
    );
    for search_idx in 0..max_match_count {
        if rc_primer_sequence.as_bytes()[search_idx]
            != forward_sequence.as_bytes()[search_start + search_idx]
        {
            return search_start + search_idx;
        }
    }
    return search_start + max_match_count;
}

fn trim_reverse_read(reverse_sequence: &str, trim_sequence: &str) -> usize {
    //
    // A match must be at least 5 bp long
    //
    let min_match_size = 5;

    //
    // Consider a seed region of 8 bp
    //
    let seed_max_length = 8;

    //
    // Allow no more than 1 mismatch in the seed region
    //
    let seed_mismatch_limit = 1;

    //
    // allow no more than 25% mismatches to entire trimmed sequence
    //
    let mismmatch_threshold_fraction = 0.25;

    if reverse_sequence.len() < min_match_size || trim_sequence.len() < min_match_size {
        return reverse_sequence.len();
    }

    'search_start: for search_start in 0..(reverse_sequence.len() - min_match_size) {
        let max_match_count = cmp::min(trim_sequence.len(), reverse_sequence.len() - search_start);
        let mismatch_limit = (max_match_count as f32 * mismmatch_threshold_fraction).floor() as u16;
        let mut mismatch_count = 0;
        for search_idx in 0..max_match_count {
            if reverse_sequence.as_bytes()[search_start + search_idx]
                != trim_sequence.as_bytes()[search_idx]
            {
                mismatch_count += 1;
            }
            if (mismatch_count > mismatch_limit)
                || (search_idx < seed_max_length && mismatch_count > seed_mismatch_limit)
            {
                continue 'search_start;
            }
        }

        //
        // There is a match. Trim aggressively without checking
        // for a better match in a latter position
        //
        return search_start;
    }

    //
    // We never found a match, don't trim anything
    //
    return reverse_sequence.len();
}

fn main() -> Result<()> {
    let args = Args::parse();

    let barcode_spacer;
    if args.solo {
        barcode_spacer = args
            .barcode_spacer
            .unwrap_or("barcode spacer must be supplied with no UMI".to_string());
    } else {
        barcode_spacer = "".to_string();
    }

    let output_dir = "output";
    let forward_paired_lines = read_lines(&args.forward_paired)?;
    let reverse_paired_lines = read_lines(&args.reverse_paired)?;

    let mut mymap: HashMap<String, Vec<String>> = HashMap::default();

    let rc_primer_sequence = revcomp(&args.pcr_primer);

    for ((header_r1, seq_r1, _sep_r1, qualities_r1), (header_r2, seq_r2, _sep_r2, qualities_r2)) in
        zip(forward_paired_lines.tuples(), reverse_paired_lines.tuples())
    {
        let seq_r1 = seq_r1?;
        let seq_r2 = seq_r2?;
        let qualities_r1 = qualities_r1?;

        if seq_r1.len() < args.min_length || seq_r2.len() < args.min_length {
            continue;
        }

        let observed_umi = &seq_r1[0..args.umi_length];

        //
        // Find the position of the (known) sample index within the forward read
        //
        let sample_index_position = match args.solo {
            true => find_solo_plate_index_position(&seq_r1, args.umi_length, &barcode_spacer), // TODO: update this to search for the position of the sample index
            false => args.umi_length,
        };

        if sample_index_position == 0 {
            continue;
        }

        //
        // the umi + sample index + primer sequence is trimmed from the 5' end of
        // the forward sequence. 'forward_trim_position' gives the amount to trim
        // from the 5' end
        //
        let forward_trim_position = trim_forward_read(
            &seq_r1,
            &rc_primer_sequence,
            sample_index_position + args.sample_index.len(),
        );

        //
        // Check if we still have enough sequence left after trimming the forward read
        //
        if seq_r1.len() - forward_trim_position < args.min_length {
            continue;
        }

        //
        // Take the reverse complement of the sequence trimmed from the 5' end of the forward read
        // and try to trim from the 3' end of reverse read
        //
        let rc_forward_trimmed_sequence = revcomp(&seq_r1[..forward_trim_position]);
        let reverse_read_trim_position = trim_reverse_read(&seq_r2, &rc_forward_trimmed_sequence);

        //
        // Check if we have enough sequence left after trimming the reverse read
        //
        if reverse_read_trim_position < args.min_length {
            continue;
        }

        let value = mymap.entry(observed_umi.to_string()).or_default();

        value.push(header_r1?);
        value.push(seq_r1[forward_trim_position..(seq_r1.len() - 2)].to_string());
        value.push(qualities_r1[forward_trim_position..(seq_r1.len() - 2)].to_string());

        value.push(header_r2?);
        value.push(seq_r2[2..(reverse_read_trim_position - 2)].to_string());
        value.push(qualities_r2?[2..(reverse_read_trim_position - 2)].to_string());
    }

    create_dir(&output_dir)?;
    for (key, value) in &mymap {
        //
        // check if there are enough reads to assemble
        //
        if value.len() > 2 * args.min_count {
            let mut handle_r1 = BufWriter::with_capacity(
                262144,
                File::create(format!("{}/{}_R1.fastq", &output_dir, &key))?,
            );
            let mut handle_r2 = BufWriter::with_capacity(
                262144,
                File::create(format!("{}/{}_R2.fastq", &output_dir, &key))?,
            );
            let mut index = 0;
            while index < value.len() {
                writeln!(handle_r1, "{}", value[index + 0])?;
                writeln!(handle_r1, "{}", value[index + 1])?;
                writeln!(handle_r1, "+")?;
                writeln!(handle_r1, "{}", value[index + 2])?;

                writeln!(handle_r2, "{}", value[index + 3])?;
                writeln!(handle_r2, "{}", value[index + 4])?;
                writeln!(handle_r2, "+")?;
                writeln!(handle_r2, "{}", value[index + 5])?;

                index += 6;
            }
        }
    }
    return Ok(());
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> anyhow::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
