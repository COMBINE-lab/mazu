use kmers::naive_impl::{CanonicalKmer, CanonicalKmerIterator};
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use super::{GetRefPos, ModIndex, U2Pos, K2U};

pub trait Validate {
    fn validate_self(&self);
    fn validate_fasta<P: AsRef<Path>>(&self, p: P);
    fn validate_ckmers(&self, ri: usize, seq: &[u8]);
}

impl<H, T> Validate for ModIndex<H, T>
where
    Self: GetRefPos,
    H: K2U,
    T: U2Pos,
{
    // Validate against internal RefSeqCollection
    fn validate_self(&self) {
        assert!(self.has_refseq());
        for r in self.iter_refs() {
            for (pos, km) in r.iter_kmers(self.k()).enumerate() {
                let km = CanonicalKmer::from(km);
                let mrps = self.get_ref_pos_eager(&km);
                if mrps.is_none() {
                    panic!(
                        "No MRPs found for true +ve kmer in ref {} @ pos {}",
                        r.ref_id(),
                        pos
                    );
                }
                let mrps = mrps.unwrap();
                let mut found = false;
                for mrp in mrps {
                    found |= (mrp.pos == pos) && (mrp.ref_id == r.ref_id());
                }

                if !found {
                    panic!(
                        "No matching MRP for true +ve kmer in ref {} @ pos {}",
                        r.ref_id(),
                        pos
                    );
                }
            }
        }
    }

    // Validate canonical k-mers of u8 slice with given reference ID
    fn validate_ckmers(&self, ref_id: usize, seq: &[u8]) {
        let k = self.k() as u8;
        let iter = CanonicalKmerIterator::from_u8_slice(seq, k);
        for tup in iter {
            let pos = tup.pos();
            let km = tup.km;
            let mrps = self.get_ref_pos_eager(&km);
            if mrps.is_none() {
                panic!(
                    "No MRPs found for true +ve kmer in ref {} @ pos {}, '{}'",
                    ref_id, pos, km
                );
            }
            let mrps = mrps.unwrap();
            let mut found = false;
            for mrp in mrps {
                found |= (mrp.pos == pos) && (mrp.ref_id == ref_id);
            }

            if !found {
                panic!(
                    "No MRPs found for true +ve kmer in ref {} @ pos {}, '{}'",
                    ref_id, pos, km
                );
            }
        }
    }

    // Validate against a fasta file
    fn validate_fasta<P: AsRef<Path>>(&self, fp: P) {
        let f = File::open(fp).unwrap();
        let f = BufReader::new(f);

        let mut lines = f.lines();

        // let mut ref_name = String::new(); // FIXME: handle ref name too?
        let mut seq_buf = Vec::new();
        let mut ref_id = -1_isize;

        loop {
            let line: Option<Result<String, std::io::Error>> = lines.next();

            if let Some(seq) = line {
                let seq = seq.unwrap();
                let seq = seq.as_bytes();
                if seq.starts_with(b">") {
                    self.validate_ckmers(ref_id as usize, &seq_buf);
                    ref_id += 1;
                    // ref_name = seq;
                    seq_buf.clear();
                } else {
                    seq_buf.extend_from_slice(seq);
                }
            } else {
                self.validate_ckmers(ref_id as usize, &seq_buf);
                break;
            }
        }
    }
}
