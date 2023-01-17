# Rust-bwa
(forked from 10XGenomics)

Rust-bwa is a Rust wrapper of the BWA api. Pass bio::io::fastq::Records in and get rust_htslib::bam::Records back.

```
extern crate bio;
extern crate bwa;

use bio::io::fastq;
use bwa::{BwaAligner, BwaReference};
use rust_htslib::bam;

fn main() {
    let aligner = BwaAligner::from_path("/Users/kg8/Downloads/rust-bwa/tests/test_ref.fa").unwrap();
    let bam_header = aligner.create_bam_header();

    let seqs = vec![
        fastq::Record::with_attrs(
            "id",
            Some("/1"),
            b"TCATTGCTTATTATGTTCATCCCGTCAACATTCAAACGGCCTGTCTCATCATGGAAGGCGCTGAATTTAC",
            b"JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ",
        ),
        fastq::Record::with_attrs(
            "id",
            Some("/2"),
            b"GGAAAACATTATTAATGGCGTCGAGCGTCCGGTTAAAGCCGCTGAATTGTTCGCGTTTACCTTGCGTGTA",
            b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        ),
    ];

    // Align the sequences. Result is nested as a Vec<Vec<bam::Record>>.
    // Each input sequence produces a Vec<bam::Record> containing its principal
    // alignment and any split or supplementary alignments.
    let alns = aligner.align_fastq_records_nested(&seqs, true).unwrap();
    
    // Write the alignments to stdout
    let mut writer = bam::Writer::from_stdout(&bam_header, bam::Format::Sam).unwrap();
    for aln in alns {
        for rec in aln {
            writer.write(&rec).unwrap();
        }
    }
}
```

Pre-built rust bindings were generated using `bindgen` for linux using the command:

```
~/.cargo/bin/bindgen \
  --no-doc-comments \
  --allowlist-function mem_align1_core \
  --allowlist-function mem_sam_pe \
  --allowlist-function mem_opt_init \
  --allowlist-function bwa_idx_load \
  --allowlist-function bwa_idx_destroy \
  --allowlist-function mem_process_seq_pe \
  --allowlist-function mem_process_seqs \
  --allowlist-function mem_align1 \
  --allowlist-function bwa_fill_scmat \
  --allowlist-var "BWA_IDX_.*" \
  --allowlist-var "MEM_F_.*" \
  wrapper.h \
  -o src/lib.rs
```

`bindgen` can be installed using `cargo install bindgen`. See the documentation [here](https://rust-lang.github.io/rust-bindgen/command-line-usage.html).
