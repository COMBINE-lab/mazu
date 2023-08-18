use kmers::naive_impl::{CanonicalKmer, CanonicalKmerIterator};
use rayon::prelude::*;
use std::{fs::File, io::BufReader, path::Path};

use crate::util::FastaReader;

use super::{GetRefPos, ModIndex, U2Pos, K2U};

pub trait Validate {
    fn validate_self(&self);
    fn validate_fasta<P: AsRef<Path>>(&self, p: P)
    where
        Self: Sync;
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

    fn validate_fasta<P: AsRef<Path>>(&self, fp: P)
    where
        Self: Sync,
    {
        let f = File::open(fp).unwrap();
        let f = BufReader::new(f);
        let rdr = FastaReader::new(f);

        // collect seqs
        let recs: Vec<(usize, Vec<u8>)> = rdr
            .into_iter()
            .map(|rec| rec.seq.as_bytes().to_vec())
            .enumerate()
            .collect();

        recs.par_iter()
            .for_each(|(ri, seq)| self.validate_ckmers(*ri, seq));
    }
}
