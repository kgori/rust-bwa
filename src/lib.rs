// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.

//! Rust-bwa provides a simple API wrapper around the BWA aligner.
//! Pass read-pair information in, and get Rust-htslib BAM records
//! back.
//!
//! ```
//! use bwa::BwaAligner;
//!
//! let bwa = BwaAligner::from_path(&"tests/test_ref.fa").unwrap();
//!
//! let r1 = b"GATGGCTGCGCAAGGGTTCTTACTGATCGCCACGTTTTTACTGGTGTTAATGGTGCTGGCGCGTCCTTTAGGCAGCGGG";
//! let q1 = b"2222222222222222222222222222222222222222222222222222222222222222222222222222222";
//! let r2 = b"TGCTGCGTAGCAGATCGACCCAGGCATTCCCTAGCGTGCTCATGCTCTGGCTGGTAAACGCACGGATGAGGGCAAAAAT";
//! let q2 = b"2222222222222222222222222222222222222222222222222222222222222222222222222222222";
//!
//! let (r1_alns, _r2_alns) = bwa.align_read_pair(b"read_name", r1, q1, r2, q2);
//! println!("r1 mapping -- tid: {}, pos: {}", r1_alns[0].tid(), r1_alns[0].pos());
//! ```
//! Align a collection of fastq Records.
//!
//! ```
//! extern crate bio;
//! extern crate rust_htslib;
//!
//! use bio::io::fastq;
//! use rust_htslib::bam;
//!
//! use bwa::BwaAligner;
//!
//! let aligner = BwaAligner::from_path(&"tests/test_ref.fa").unwrap();
//! let bam_header = aligner.create_bam_header();
//!
//! let seqs = vec![
//!     fastq::Record::with_attrs(
//!         "id",
//!         Some("/1"), // optional description, ignored by BWA
//!         b"TCATTGCTTATTATGTTCATCCCGTCAACATTCAAACGGCCTGTCTCATCATGGAAGGCGCTGAATTTAC",
//!         b"JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ",
//!     ),
//!     fastq::Record::with_attrs(
//!         "id",
//!         Some("/2"),
//!         b"GGAAAACATTATTAATGGCGTCGAGCGTCCGGTTAAAGCCGCTGAATTGTTCGCGTTTACCTTGCGTGTA",
//!         b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
//!         ),
//!     ];
//!
//! let alns = aligner
//!     .align_fastq_records_nested(&seqs, /* paired */ true)
//!     .unwrap();
//! assert!(alns[0][0].pos() == 630);
//! assert!(alns[1][0].pos() == 700);
//! ```

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

extern crate bio;
extern crate libc;
extern crate rust_htslib;

extern crate bwa_sys;
extern crate thiserror;

use std::ffi::{c_char, CStr, CString};
use std::path::Path;
use std::sync::{Arc, Mutex};

use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::record::Record;
use rust_htslib::bam::HeaderView;

use bio::io::fastq;
use bwa_sys::bseq1_t;

// include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

/// BWA settings object. Currently only default settings are enabled
pub struct BwaSettings {
    bwa_settings: bwa_sys::mem_opt_t,
}

impl BwaSettings {
    /// Create a `BwaSettings` object with default BWA parameters
    pub fn new() -> BwaSettings {
        let ptr = unsafe { bwa_sys::mem_opt_init() };
        let bwa_settings = unsafe { *ptr };
        unsafe { libc::free(ptr as *mut libc::c_void) };
        BwaSettings { bwa_settings }
    }

    /// Set alignment scores
    pub fn set_scores(
        mut self,
        matchp: i32,
        mismatch: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> BwaSettings {
        self.bwa_settings.a = matchp;
        self.bwa_settings.b = mismatch;
        self.bwa_settings.o_del = gap_open;
        self.bwa_settings.o_ins = gap_open;
        self.bwa_settings.e_del = gap_extend;
        self.bwa_settings.e_ins = gap_extend;

        unsafe {
            bwa_sys::bwa_fill_scmat(matchp, mismatch, self.bwa_settings.mat.as_mut_ptr());
        }
        self
    }

    /// Set clipping score penalties
    pub fn set_clip_scores(mut self, clip5: i32, clip3: i32) -> BwaSettings {
        self.bwa_settings.pen_clip5 = clip5;
        self.bwa_settings.pen_clip3 = clip3;
        self
    }

    /// Set unpaired read penalty
    pub fn set_unpaired(mut self, unpaired: i32) -> BwaSettings {
        self.bwa_settings.pen_unpaired = unpaired;
        self
    }

    /// Mark shorter splits as secondary
    pub fn set_no_multi(mut self) -> BwaSettings {
        self.bwa_settings.flag |= 0x10; // MEM_F_NO_MULTI
        self
    }
}

#[derive(Debug, thiserror::Error)]
#[error("{0}")]
pub struct ReferenceError(String);

#[derive(Debug, thiserror::Error)]
#[error("{0}")]
pub struct BwaAlignmentError(String);

/// A BWA reference object to perform alignments to.
/// Must be loaded from a BWA index created with `bwa index`
pub struct BwaReference {
    bwt_data: *const bwa_sys::bwaidx_t,
    contig_names: Vec<String>,
    contig_lengths: Vec<usize>,
}
unsafe impl Sync for BwaReference {}

impl BwaReference {
    /// Load a BWA reference from disk. Pass the fasta filename of the
    /// original reference as `path`
    pub fn open<P: AsRef<Path>>(path: P) -> Result<BwaReference, ReferenceError> {
        let idx_file = CString::new(path.as_ref().to_str().unwrap()).unwrap();
        let idx = unsafe { bwa_sys::bwa_idx_load(idx_file.as_ptr(), 0x7 as i32) }; // FIXME -- use BWA_IDX_ALL

        if idx.is_null() {
            return Err(ReferenceError(format!(
                "couldn't load reference: {:?}",
                path.as_ref()
            )));
        }

        let mut contig_names = Vec::new();
        let mut contig_lengths = Vec::new();
        let num_contigs = unsafe { (*(*idx).bns).n_seqs };

        for i in 0..num_contigs as isize {
            unsafe {
                let name = CStr::from_ptr((*(*(*idx).bns).anns.offset(i)).name);
                let sz = (*(*(*idx).bns).anns.offset(i)).len;

                let name_string = name.to_owned().into_string().unwrap();
                contig_names.push(name_string);
                contig_lengths.push(sz as usize)
            }
        }

        Ok(BwaReference {
            bwt_data: idx,
            contig_names,
            contig_lengths,
        })
    }

    pub fn create_bam_header(&self) -> Header {
        let mut header = Header::new();
        self.populate_bam_header(&mut header);
        header
    }

    pub fn populate_bam_header(&self, header: &mut Header) {
        for (ref contig_name, &len) in self.contig_names.iter().zip(self.contig_lengths.iter()) {
            add_ref_to_bam_header(header, &contig_name, len);
        }
    }
}

impl Drop for BwaReference {
    fn drop(&mut self) {
        unsafe {
            bwa_sys::bwa_idx_destroy(self.bwt_data as *mut bwa_sys::bwaidx_t);
        }
    }
}

fn add_ref_to_bam_header(header: &mut Header, seq_name: &str, seq_len: usize) {
    let mut header_rec = HeaderRecord::new(b"SQ");
    header_rec.push_tag(b"SN", &seq_name);
    header_rec.push_tag(b"LN", &seq_len);
    header.push_record(&header_rec);
}

/// Paired-end statistics structure used by BWA to score paired-end reads
pub struct PairedEndStats {
    inner: [bwa_sys::mem_pestat_t; 4],
}

impl PairedEndStats {
    /// Generate a 'simple' paired-end read structure that standard forward-reverse
    /// pairs as created by TruSeq, Nextera, or Chromium Genome sample preparations.
    pub fn simple(avg: f64, std: f64, low: i32, high: i32) -> PairedEndStats {
        let pe_stat_null = || bwa_sys::mem_pestat_t {
            failed: 1,
            low: 0,
            high: 0,
            avg: 0.0,
            std: 100.0,
        };

        let pes = [
            pe_stat_null(),
            bwa_sys::mem_pestat_t {
                failed: 0,
                low,
                high,
                avg,
                std,
            },
            pe_stat_null(),
            pe_stat_null(),
        ];

        PairedEndStats { inner: pes }
    }

    pub fn default() -> PairedEndStats {
        Self::simple(200.0, 100.0, 35, 600)
    }
}

/// A BWA aligner. Carries everything required to align
/// reads to a reference and generate BAM records.
pub struct BwaAligner {
    reference: BwaReference,
    header_view: Arc<Mutex<HeaderView>>,
    settings: BwaSettings,
    pe_stats: PairedEndStats,
}
// this is not automatically derived because of an interior
//   mutable pointer inside HeaderView. It _is_ mutated
//   by the Record::from_sam function, so guard it with a mutex
unsafe impl Sync for BwaAligner {}

impl BwaAligner {
    /// Load a BWA reference from the given path and use default BWA settings and paired-end structure.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<BwaAligner, ReferenceError> {
        let bwa_ref = BwaReference::open(path)?;
        Ok(BwaAligner::new(
            bwa_ref,
            BwaSettings::new(),
            PairedEndStats::default(),
        ))
    }

    pub fn new(
        reference: BwaReference,
        settings: BwaSettings,
        pe_stats: PairedEndStats,
    ) -> BwaAligner {
        let header = reference.create_bam_header();
        let header_view = Arc::new(Mutex::new(HeaderView::from_header(&header)));

        BwaAligner {
            reference,
            header_view,
            settings,
            pe_stats,
        }
    }

    pub fn create_bam_header(&self) -> Header {
        self.reference.create_bam_header()
    }

    /// Align an array of `bio::io::fastq::Records` to the reference and return a `Vec<bool>`
    /// to indicate whether each `fastq::Record` mapped to the reference (true) or was
    /// unmapped (false).
    pub fn get_alignment_status(
        &self,
        records: &[fastq::Record],
        paired: bool,
    ) -> Result<Vec<bool>, BwaAlignmentError> {
        let (cnames, mut vseqs, mut vquals) = Self::extract_fastqs(records);
        let mut bseqs = BseqVec::new(cnames, &mut vseqs, &mut vquals);
        self.align_bseqs(paired, &mut bseqs);

        if let Some(sams) = bseqs.to_sams() {
            Ok(sams
                .iter()
                .map(|s| {
                    self.parse_sam_to_records(s.to_bytes())
                        .iter()
                        .any(|r| !r.is_unmapped())
                })
                .collect())
        } else {
            Err(BwaAlignmentError("An alignment error occurred".to_string()))
        }
    }

    fn validate_paired_records(records: &[fastq::Record], paired: bool) -> Result<(), BwaAlignmentError> {
        if paired {
            if records.len() & 1 == 1 {
                return Err(BwaAlignmentError("Expected an even number of paired reads".to_string()));
            }
            if records.chunks(2).any(|a| a[0].id() != a[1].id()) {
                return Err(BwaAlignmentError("Paired read names don't match".to_string()));
            }
        }
        Ok(())
    }

    /// Align an array of `bio::io::fastq::Records` to the reference and return a Vec
    /// of alignments as `rust_htslib::bam::Records`.
    pub fn align_fastq_records(
        &self,
        records: &[fastq::Record],
        paired: bool,
    ) -> Result<Vec<Record>, BwaAlignmentError> {
        Self::validate_paired_records(records, paired)?;
        let (cnames, mut vseqs, mut vquals) = Self::extract_fastqs(records);
        let mut bseqs = BseqVec::new(cnames, &mut vseqs, &mut vquals);
        self.align_bseqs(paired, &mut bseqs);
        if let Some(sams) = bseqs.to_sams() {
            Ok(self.sams_to_bam_records(sams))
        } else {
            Err(BwaAlignmentError("An alignment error occurred".to_string()))
        }
    }

    /// Align an array of `bio::io::fastq::Records` records to the reference and return a Vec
    /// of subalignments for each `fastq::record` as a Vec of `bam::Records`.
    pub fn align_fastq_records_nested(
        &self,
        records: &[fastq::Record],
        paired: bool,
    ) -> Result<Vec<Vec<Record>>, BwaAlignmentError> {
        Self::validate_paired_records(records, paired)?;
        let (cnames, mut vseqs, mut vquals) = Self::extract_fastqs(records);
        let mut bseqs = BseqVec::new(cnames, &mut vseqs, &mut vquals);
        self.align_bseqs(paired, &mut bseqs);
        if let Some(sams) = bseqs.to_sams() {
            Ok(sams
                .iter()
                .map(|s| self.parse_sam_to_records(s.to_bytes()))
                .collect())
        } else {
            Err(BwaAlignmentError("An alignment error occurred".to_string()))
        }
    }

    fn sams_to_bam_records(&self, sams: Vec<&CStr>) -> Vec<Record> {
        sams.iter()
            .flat_map(|s| self.parse_sam_to_records(s.to_bytes()))
            .collect()
    }

    fn align_bseqs(&self, paired: bool, bseqs: &mut BseqVec) {
        unsafe {
            let r = *(self.reference.bwt_data);
            let mut settings = self.settings.bwa_settings;
            if paired {
                settings.flag |= bwa_sys::MEM_F_PE as i32;
            } else {
                settings.flag &= !(bwa_sys::MEM_F_PE as i32);
            }
            bwa_sys::mem_process_seqs(
                &settings,
                r.bwt,
                r.bns,
                r.pac,
                0,
                bseqs.inner.len() as i32,
                bseqs.inner.as_mut_ptr(),
                self.pe_stats.inner.as_ptr(),
            );
        }
        bseqs.aligned = true;
    }

    fn extract_fastqs(records: &[fastq::Record]) -> (Vec<*mut c_char>, Vec<Vec<u8>>, Vec<Vec<u8>>) {
        let cnames = records
            .iter()
            .map(|rec| CString::new(rec.id()).unwrap().into_raw())
            .collect::<Vec<_>>();
        let vseqs = records
            .iter()
            .map(|rec| Vec::from(rec.seq()))
            .collect::<Vec<_>>();
        let vquals = records
            .iter()
            .map(|rec| Vec::from(rec.qual()))
            .collect::<Vec<_>>();
        (cnames, vseqs, vquals)
    }

    /// Align a read-pair to the reference.
    pub fn align_read_pair(
        &self,
        name: &[u8],
        r1: &[u8],
        q1: &[u8],
        r2: &[u8],
        q2: &[u8],
    ) -> (Vec<Record>, Vec<Record>) {
        let read1 = fastq::Record::with_attrs(std::str::from_utf8(name).unwrap(), None, r1, q1);

        let read2 = fastq::Record::with_attrs(std::str::from_utf8(name).unwrap(), None, r2, q2);

        let reads = vec![read1, read2];
        let aln = self.align_fastq_records_nested(&reads, true).unwrap();
        let mut it = aln.iter();
        (it.next().unwrap().to_vec(), it.next().unwrap().to_vec())
    }

    fn parse_sam_to_records(&self, sam: &[u8]) -> Vec<Record> {
        let mut records = Vec::new();

        for slc in sam.split(|x| *x == b'\n') {
            if slc.len() > 0 {
                let record = {
                    let header_view = self.header_view.lock().unwrap();
                    Record::from_sam(&header_view, slc).unwrap()
                };
                records.push(record);
            }
        }

        records
    }
}

/// A wrapper around a vector of the BWA sequence type bseq1_t
struct BseqVec {
    inner: Vec<bseq1_t>,
    aligned: bool,
}

impl BseqVec {
    fn new(names: Vec<*mut c_char>, seqs: &mut [Vec<u8>], quals: &mut [Vec<u8>]) -> BseqVec {
        let mut bseqs = Vec::new();
        for i in 0..names.len() {
            let bseq = bwa_sys::bseq1_t {
                l_seq: seqs[i].len() as i32,
                name: names[i],
                seq: seqs[i].as_mut_ptr() as *mut i8,
                qual: quals[i].as_mut_ptr() as *mut i8,
                comment: std::ptr::null_mut(),
                id: i as i32,
                sam: std::ptr::null_mut(),
            };
            bseqs.push(bseq);
        }
        BseqVec {
            inner: bseqs,
            aligned: false,
        }
    }

    /// Extracts the SAM information from the aligned bseq1_ts as a vec of &CStr
    fn to_sams(&self) -> Option<Vec<&CStr>> {
        if !self.aligned {
            return None;
        }
        let sams = self
            .inner
            .iter()
            .map(|b| unsafe { CStr::from_ptr(b.sam) })
            .collect::<Vec<_>>();
        Some(sams)
    }
}

impl Drop for BseqVec {
    fn drop(&mut self) {
        for bseq in &self.inner {
            unsafe {
                libc::free(bseq.name as *mut libc::c_void);
                libc::free(bseq.sam as *mut libc::c_void);
                libc::free(bseq.comment as *mut libc::c_void);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    extern crate itertools;
    use tests::itertools::Itertools;

    fn load_aligner() -> BwaAligner {
        let aln = BwaAligner::from_path("tests/test_ref.fa");
        aln.unwrap()
    }

    #[test]
    fn test_load_aligner() {
        let _ = load_aligner();
    }

    fn read_simple() -> [&'static [u8]; 5] {
        let name: &[u8] = b"@chr_727436_727956_3:0:0_1:0:0_0/1";
        let r1  : &[u8] = b"GATGGCTGCGCAAGGGTTCTTACTGATCGCCACGTTTTTACTGGTGTTAATGGTGCTGGCGCGTCCTTTAGGCAGCGGGCTGGCGCGGCTGATTAATGACATTCCTCTTCCCGGTACAACGGGCGTTGAGCGCGAACTTTTTCGCGCACT";
        let q1  : &[u8] = b"222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222";

        let r2  : &[u8] = b"TGCTGCGTAGCAGATCGACCCAGGCATTCCCTAGCGTGCTCATGCTCTGGCTGGTAAACGCACGGATGAGGGCAAAAATCACCGCAATCCCGCTGGCGGCAGAAAGAAAGTTTTGCACCGTTAAGCCCGCCATCTGGCTGAAATAGCTCA";
        let q2  : &[u8] = b"222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222";

        [name, r1, q1, r2, q2]
    }

    fn read_split() -> [&'static [u8]; 5] {
        let name = b"@chr_1561275_1561756_1:0:0_2:0:0_5c/1";
        let r1 = b"GCATCGATAAGCAGGTCAAATTCTCCCGTCATTATCACCTCTGCTACTTAAATTTCCCGCTTTATAAGCCGATTACGGCCTGGCATTACCCTATCCATAATTTAGGTGGGATGCCCGGTGCGTGGTTGGCAGATCCGCTGTTCTTTATTT";
        let q1 = b"222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222";

        let r2 = b"TCATCGACCCAGGTATCATCGCGACGGGTACGATTACTGGCGAAGGTGAGAATGTTTAAAATCCAGCCGCCGAGTTTTTCAGCAATGGTCACCCATGACCAACCGGTGAACAACGTGAGGGCCGCTGCCCAAACGCATAGCAGCGCAATA";
        let q2 = b"222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222";

        [name, r1, q1, r2, q2]
    }

    fn align_read(r: [&[u8]; 5]) -> (Vec<Record>, Vec<Record>) {
        let bwa = load_aligner();
        bwa.align_read_pair(r[0], r[1], r[2], r[3], r[4])
    }

    #[test]
    fn simple_align() {
        let (r1, r2) = align_read(read_simple());
        assert_eq!(r1[0].pos(), 727806);
        assert_eq!(r2[0].pos(), 727435);
    }

    #[test]
    fn split_align() {
        let (r1, r2) = align_read(read_split());
        assert_eq!(r1.len(), 2);
        assert_eq!(r1[0].pos(), 931375);
        assert_eq!(r1[1].pos(), 932605);
        assert_eq!(r2[0].pos(), 932937);
    }

    #[test]
    fn header() {
        let reference = BwaReference::open("tests/test_ref.fa").unwrap();
        let hdr = b"@SQ\tSN:PhiX\tLN:5386\n@SQ\tSN:chr\tLN:4639675";
        assert_eq!(
            reference.create_bam_header().to_bytes().as_slice(),
            &hdr[..]
        );
    }

    fn read_fastq() -> Result<Vec<fastq::Record>, bio::io::fastq::Error> {
        let reader1 = fastq::Reader::from_file("tests/test.1.fq").unwrap();
        let reader2 = fastq::Reader::from_file("tests/test.2.fq").unwrap();
        let records: Result<Vec<fastq::Record>, _> =
            reader1.records().interleave(reader2.records()).collect();
        records
    }

    #[test]
    fn test_align_fastqs() {
        let records = read_fastq().unwrap();
        let bwa = load_aligner();
        let recs = bwa.align_fastq_records(&records, true).unwrap();
        assert_eq!(recs[0].pos(), 1330);
    }

    #[test]
    fn test_do_reads_align() {
        let records = read_fastq().unwrap();
        let bwa = load_aligner();
        let result = bwa.get_alignment_status(&records, true).unwrap();
        assert!(result[0]);
        assert!(result[1]);
        assert!(result[2]);
        assert!(result[3]);
        assert!(result[4]);
        assert!(result[5]);
        assert!(result[6]);
        assert!(result[7]);
        assert!(result[8]);
        assert!(result[9]);
        assert!(result[10]);
        assert!(result[11]);
        assert!(!result[12]);
        assert!(!result[13]);
    }

    #[test]
    fn test_odd_paired_reads_throw_error() {
        let records = &read_fastq().unwrap()[..3];
        let bwa = load_aligner();
        assert!(bwa.align_fastq_records(&records, true).is_err());
    }

    #[test]
    fn test_unordered_reads_throw_error() {
        let mut records_1 = read_fastq().unwrap().into_iter().step_by(2).collect::<Vec<fastq::Record>>();
        let mut records_2 = read_fastq().unwrap().into_iter().skip(1).step_by(2).collect::<Vec<fastq::Record>>();
        records_1.append(&mut records_2);
        let bwa = load_aligner();
        assert!(bwa.align_fastq_records(records_1.as_slice(), true).is_err());
    }

    #[test]
    fn test_aligner_can_generate_header() {
        let bwa = load_aligner();
        let hdr = b"@SQ\tSN:PhiX\tLN:5386\n@SQ\tSN:chr\tLN:4639675";
        assert_eq!(
            bwa.create_bam_header().to_bytes().as_slice(),
            &hdr[..]
        );
    }
}
