use anyhow::{Ok, Result};
use clap::Parser;
use core::str::from_utf8;
use itertools::Itertools;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::fs;
use std::fs::{create_dir, File};
use std::io::{self, BufRead, BufWriter, Write};
use std::iter::zip;
use std::path::Path;

use flate2::write::GzEncoder;
use flate2::Compression;

#[derive(Parser, Debug)]
#[clap(author="Ryan Kelley", version, about="Demultiplex loop samples", long_about = None)]
struct Args {
    /// Uncompressed FASTQ file with paired forward reads
    forward_paired: String,

    /// Uncompressed FASTQ file with paired reverse reads
    reverse_paired: String,

    /// Uncompressed FASTQ file with unpaired forward reads
    forward_unpaired: String,

    /// comma-delimited list of barcodes for de-multiplexing
    barcodes: String,

    /// Nucleotide sequence of PCR primer
    #[clap(long, default_value = "")]
    pcr_primer: String,

    /// Amount of 5' context to use from PCR primer for barcode mapping
    #[clap(long, default_value_t = 2)]
    pcr_primer_context_length: usize,

    /// Length of UMI prefix to use for UMI grouping
    #[clap(long, default_value_t = 2)]
    umi_group_length: usize,

    /// Length of the UMI
    #[clap(long, default_value_t = 16)]
    umi_length: usize,

    /// Enable processing of Solo data
    #[clap(long, takes_value = false)]
    solo: bool,

    /// Nucleotide sequence of barcode spacer
    #[clap(long)]
    barcode_spacer: Option<String>,
}

/// Clean up empty files
fn clean_files(
    indices: &Vec<String>,
    pcr_primer_context_length: usize,
    umi_group_length: usize,
) -> Result<()> {
    for index in indices.iter() {
        let directory = format!(
            "sample_{}",
            &index[..(index.len() - pcr_primer_context_length)]
        );
        for umi_bin in generate_all_strings("ACGT", umi_group_length) {
            let r1_filename = format!(
                "{}/{}_{}_R1.fastq.gz",
                &directory,
                &index[..(index.len() - pcr_primer_context_length)],
                &umi_bin
            );
            let r2_filename = format!(
                "{}/{}_{}_R2.fastq.gz",
                &directory,
                &index[..(index.len() - pcr_primer_context_length)],
                &umi_bin
            );
            let unpaired_filename = format!("{}/{}_R1_unpaired.fastq.gz", &directory, &umi_bin);
            let metadata = fs::metadata(r1_filename.to_string());
            let metadata = metadata?;
            println!("size of {} is {}", &r1_filename, metadata.len());
            if metadata.len() <= 4096 {
                println!("clean up {}", &r1_filename);
                let _ = std::fs::remove_file(&r1_filename);
                let _ = std::fs::remove_file(r2_filename);
                let _ = std::fs::remove_file(unpaired_filename);
            }
        }
    }
    return Ok(());
}

/// Generate compressed file handles for paired R1 and R2 output from a set of indices
fn generate_paired_file_handles_for_indices(
    indices: &Vec<String>,
    pcr_primer_context_length: usize,
    umi_group_length: usize,
) -> Result<HashMap<String, (Box<dyn Write>, Box<dyn Write>)>> {
    let mut result: HashMap<String, (Box<dyn Write>, Box<dyn Write>)> = HashMap::default();
    for index in indices.iter() {
        let directory = format!(
            "sample_{}",
            &index[..(index.len() - pcr_primer_context_length)]
        );
        create_dir(&directory)?;
        for umi_bin in generate_all_strings("ACGT", umi_group_length) {
            result.insert(
                index.to_string() + &umi_bin,
                (
                    Box::new(GzEncoder::new(
                        File::create(format!(
                            "{}/{}_{}_R1.fastq.gz",
                            &directory,
                            &index[..(index.len() - pcr_primer_context_length)],
                            &umi_bin
                        ))?,
                        Compression::default(),
                    )),
                    Box::new(GzEncoder::new(
                        File::create(format!(
                            "{}/{}_{}_R2.fastq.gz",
                            &directory,
                            &index[..(index.len() - pcr_primer_context_length)],
                            &umi_bin
                        ))?,
                        Compression::default(),
                    )),
                ),
            );
        }
    }
    return Ok(result);
}

/// Generate compressed file handles for unpaired R1 output from a set of indices
fn generate_unpaired_file_handles_for_indices(
    indices: &Vec<String>,
    pcr_primer_context_length: usize,
    umi_group_length: usize,
) -> Result<HashMap<String, Box<dyn Write>>> {
    let mut result: HashMap<String, Box<dyn Write>> = HashMap::default();
    for barcode in indices.iter() {
        let directory = format!(
            "sample_{}",
            &barcode[..(barcode.len() - pcr_primer_context_length)]
        );
        for umi_bin in generate_all_strings("ACGT", umi_group_length) {
            result.insert(
                barcode.to_string() + &umi_bin,
                Box::new(GzEncoder::new(
                    File::create(format!(
                        "{}/{}_{}_R1_unpaired.fastq.gz",
                        &directory,
                        &barcode[..(barcode.len() - pcr_primer_context_length)],
                        &umi_bin
                    ))?,
                    Compression::default(),
                )),
            );
        }
    }
    return Ok(result);
}

fn find_solo_plate_index<'a, 'b>(
    seq: &'a str,
    bc_len: usize,
    spacer: &'b str,
    index_len: usize,
) -> &'a str {
    for idx in 0..4 {
        let spacer_start = bc_len + idx + 1;
        let spacer_end = spacer_start + spacer.len();
        let index_end = spacer_start + spacer.len() + index_len;
        if index_end < seq.len() && &seq[spacer_start..spacer_end] == spacer {
            return &seq[spacer_end..index_end];
        }
    }
    return &"";
}

fn generate_all_strings(alphabet: &str, n: usize) -> Vec<String> {
    if n == 0 {
        return vec!["".to_string()];
    }

    let mut strings = Vec::new();
    for char in alphabet.chars() {
        for string in generate_all_strings(alphabet, n - 1) {
            strings.push(char.to_string() + &string);
        }
    }
    return strings;
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

fn revcomp(dna: &str) -> String {
    // result vector
    let mut rdna: String = String::with_capacity(dna.len());

    // iterate through the input &str
    for c in dna.chars().rev() {
        rdna.push(switch_base(c));
    }
    return rdna;
}

fn main() -> Result<()> {
    let args = Args::parse();

    let rc_primer = revcomp(&args.pcr_primer);

    //
    // Get information about indices
    //
    let original_indices = args.barcodes.split(",").collect_vec();

    //
    // Add 5' context from PCR primer into the sample index sequences
    //
    let mut indices: Vec<String> = Vec::new();
    for index in original_indices {
        let index_with_primer_context = &[
            index.as_bytes(),
            &rc_primer.as_bytes()[..args.pcr_primer_context_length],
        ]
        .concat();
        indices.push(from_utf8(index_with_primer_context).unwrap().to_string());
    }

    {
        let first_index_length = indices
            .iter()
            .next()
            .expect("List of provided indices is empty")
            .len();
        //
        // Construct output file handles
        //
        let mut paired_file_handles = generate_paired_file_handles_for_indices(
            &indices,
            args.pcr_primer_context_length,
            args.umi_group_length,
        )?;
        let mut unpaired_file_handles = generate_unpaired_file_handles_for_indices(
            &indices,
            args.pcr_primer_context_length,
            args.umi_group_length,
        )?;

        //
        // Construct counts dictionary
        //
        let mut index_counts: HashMap<String, u32> = HashMap::default();
        for index in indices.iter() {
            index_counts.entry(index.to_string()).or_insert(0);
        }

        let forward_paired_lines = read_lines(&args.forward_paired)?;
        let reverse_paired_lines = read_lines(&args.reverse_paired)?;

        let barcode_spacer;
        if args.solo {
            barcode_spacer = args
                .barcode_spacer
                .unwrap_or("barcode spacer must be supplied with no UMI".to_string());
        } else {
            barcode_spacer = "".to_string();
        }

        for (
            (header_r1, seq_r1, sep_r1, qualities_r1),
            (header_r2, seq_r2, sep_r2, qualities_r2),
        ) in zip(forward_paired_lines.tuples(), reverse_paired_lines.tuples())
        {
            let seq_r1 = seq_r1?;

            if seq_r1.len() > args.umi_length + first_index_length {
                let observed_barcode;
                if args.solo {
                    observed_barcode = find_solo_plate_index(
                        &seq_r1,
                        args.umi_length,
                        &barcode_spacer,
                        first_index_length,
                    );
                } else {
                    observed_barcode =
                        &seq_r1[args.umi_length..(args.umi_length + first_index_length)];
                }
                let file_handles_for_barcode = paired_file_handles
                    .get_mut(&(observed_barcode.to_owned() + &seq_r1[0..args.umi_group_length]));
                if let Some((r1_handle, r2_handle)) = file_handles_for_barcode {
                    writeln!(r1_handle, "{}", header_r1?)?;
                    writeln!(r1_handle, "{}", seq_r1)?;
                    writeln!(r1_handle, "{}", sep_r1?)?;
                    writeln!(r1_handle, "{}", qualities_r1?)?;

                    writeln!(r2_handle, "{}", header_r2?)?;
                    writeln!(r2_handle, "{}", seq_r2?)?;
                    writeln!(r2_handle, "{}", sep_r2?)?;
                    writeln!(r2_handle, "{}", qualities_r2?)?;
                }

                let index_counts_length = index_counts.len();
                match index_counts.entry(observed_barcode.to_owned()) {
                    Entry::Occupied(mut o) => {
                        o.insert(o.get() + 1);
                    }
                    Entry::Vacant(v) => {
                        if index_counts_length < 1000 {
                            v.insert(1);
                        }
                    }
                }
            }
        }

        let forward_unpaired_lines = read_lines(&args.forward_unpaired)?;
        for (header_r1, seq_r1, sep_r1, qualities_r1) in forward_unpaired_lines.tuples() {
            let seq_r1 = seq_r1?;
            if seq_r1.len() > args.umi_length + first_index_length {
                let observed_barcode;
                if args.solo {
                    observed_barcode = find_solo_plate_index(
                        &seq_r1,
                        args.umi_length,
                        &barcode_spacer,
                        first_index_length,
                    );
                } else {
                    observed_barcode =
                        &seq_r1[args.umi_length..(args.umi_length + first_index_length)];
                }
                let file_handle_for_barcode = unpaired_file_handles
                    .get_mut(&(observed_barcode.to_owned() + &seq_r1[0..args.umi_group_length]));
                if let Some(handle_r1) = file_handle_for_barcode {
                    writeln!(handle_r1, "{}", header_r1?)?;
                    writeln!(handle_r1, "{}", seq_r1)?;
                    writeln!(handle_r1, "{}", sep_r1?)?;
                    writeln!(handle_r1, "{}", qualities_r1?)?;
                }

                let index_counts_length = index_counts.len();
                match index_counts.entry(observed_barcode.to_owned()) {
                    Entry::Occupied(mut o) => {
                        o.insert(o.get() + 1);
                    }
                    Entry::Vacant(v) => {
                        if index_counts_length < 1000 {
                            v.insert(1);
                        }
                    }
                }
            }
        }

        //
        // Write out the demult_summary.txt
        //
        let mut index_counts_writer =
            BufWriter::with_capacity(262144, File::create("demult_summary.txt")?);
        writeln!(
            index_counts_writer,
            "sample_id\tindex_seq\tdemult_read_count"
        )?;
        for (a, b) in indices.iter().enumerate() {
            writeln!(
                index_counts_writer,
                "sample_{}\t{}\t{}",
                a + 1,
                b,
                index_counts[b]
            )?;
        }
    }
    clean_files(
        &indices,
        args.pcr_primer_context_length,
        args.umi_group_length,
    )?;
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
