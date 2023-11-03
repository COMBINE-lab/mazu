use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use kmers::naive_impl::{
    seq_vector::{SeqVecKmerIterator, SeqVector, SeqVectorSlice},
    Base, Kmer,
};
use serde::{Deserialize, Serialize};

use crate::spt::SPT;

// FIXME: need to fix to account for polyNs for contig sampling where unitig occurrences are iterated over
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RefSeqCollection {
    pub(crate) seq: Option<SeqVector>,
    pub(crate) prefix_sum: Vec<usize>, // prefix sum of reference sequence lengths
}

impl RefSeqCollection {
    pub fn dummy() -> Self {
        log::warn!("FIXME: USING DUMMY REF SEQ");
        Self {
            seq: None,
            prefix_sum: Vec::new(),
        }
    }

    pub fn from_spt(spt: &SPT) -> Self {
        Self {
            seq: None,
            prefix_sum: crate::util::prefix_sum(&spt.ref_lens),
        }
    }

    pub fn from_fasta<P: AsRef<Path>>(p: P) -> std::io::Result<Self> {
        // FIXME: use util::FastaReader instead
        log::warn!("Note that 'A' is inserted for each 'N' in any polyN stretch");
        let f = File::open(p)?;
        let f = BufReader::new(f);
        let mut buf = String::new();
        let mut seq = SeqVector::new();
        let mut lines = f.lines();

        let first = lines.next().expect("Empty Fasta")?;
        assert!(first.starts_with('>'));

        let mut accum = 0;
        let mut prefix_sum = vec![0];

        for line in lines {
            let mut line = line?;

            if line.starts_with('>') {
                // is sequence
                // flush and do work
                accum += buf.len();
                prefix_sum.push(accum);
                seq.push_chars(buf.as_bytes());
                buf.clear();
            } else {
                line = line.to_ascii_uppercase();
                line = line.replace('N', "A");
                buf += &line;
            }
        }

        accum += buf.len();
        prefix_sum.push(accum);
        seq.push_chars(buf.as_bytes());
        buf.clear();

        let seq = Some(seq);

        Ok(Self { seq, prefix_sum })
    }

    pub fn from_fasta_w_min_len<P: AsRef<Path>>(p: P, min_len: usize) -> std::io::Result<Self> {
        // FIXME: use util::FastaReader instead
        log::warn!("Note that 'A' is inserted for each 'N' in any polyN stretch");
        let f = File::open(p)?;
        let f = BufReader::new(f);
        let mut buf = String::new();
        let mut seq = SeqVector::new();
        let mut lines = f.lines();

        let first = lines.next().expect("Empty Fasta")?;
        assert!(first.starts_with('>'));

        let mut accum = 0;
        let mut prefix_sum = vec![0];

        for line in lines {
            let mut line = line?;

            if line.starts_with('>') {
                // is sequence
                if buf.len() >= min_len {
                    // record ref seq in buffer
                    accum += buf.len();
                    prefix_sum.push(accum);
                    seq.push_chars(buf.as_bytes());
                }
                buf.clear(); // flush
            } else {
                line = line.to_ascii_uppercase();
                line = line.replace('N', "A");
                buf += &line;
            }
        }

        if buf.len() >= min_len {
            // record ref seq in buffer
            accum += buf.len();
            prefix_sum.push(accum);
            seq.push_chars(buf.as_bytes());
        }
        buf.clear(); // flush
        let seq = Some(seq);

        Ok(Self { seq, prefix_sum })
    }
}

impl RefSeqCollection {
    pub fn from_parts(seq: Option<SeqVector>, prefix_sum: Vec<usize>) -> Self {
        Self { seq, prefix_sum }
    }

    pub fn has_seq(&self) -> bool {
        self.seq.is_some()
    }

    pub fn total_len(&self) -> usize {
        self.prefix_sum[self.n_refs()]
    }

    pub fn n_refs(&self) -> usize {
        self.prefix_sum.len() - 1
    }

    pub fn ref_len(&self, ref_id: usize) -> usize {
        assert!(ref_id < self.n_refs());
        self.prefix_sum[ref_id + 1] - self.prefix_sum[ref_id]
    }

    pub fn get_refseq(&self, i: usize) -> RefSeqSlice<'_> {
        if let Some(seq) = &self.seq {
            let s = self.prefix_sum[i];
            let e = self.prefix_sum[i + 1];

            let seq = seq.slice(s, e);
            RefSeqSlice { ref_id: i, seq }
        } else {
            panic!("Refseq is None")
        }
    }

    pub fn iter_refs(&self) -> RefSeqIter {
        RefSeqIter {
            ref_id: 0,
            refs: self,
        }
    }
}

#[derive(Debug, Clone)]
pub struct RefSeqSlice<'a> {
    ref_id: usize,
    seq: SeqVectorSlice<'a>,
}

impl<'a> ToString for RefSeqSlice<'a> {
    fn to_string(&self) -> String {
        self.seq.to_string()
    }
}

impl RefSeqSlice<'_> {
    pub fn ref_id(&self) -> usize {
        self.ref_id
    }

    pub fn iter_kmers(&self, k: usize) -> SeqVecKmerIterator {
        self.seq.iter_kmers(k)
    }

    pub fn get_kmer(&self, k: usize, pos: usize) -> Kmer {
        self.seq.get_kmer(pos, k)
    }

    pub fn get_base(&self, pos: usize) -> Base {
        self.seq.get_base(pos)
    }

    pub fn get_base_checked(&self, pos: usize) -> Option<Base> {
        if pos < self.len() {
            Some(self.get_base(pos))
        } else {
            None
        }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.len() == 0
    }
}

// Iterator over references
#[derive(Debug)]
pub struct RefSeqIter<'a> {
    ref_id: usize,
    refs: &'a RefSeqCollection,
}

impl<'a> RefSeqIter<'a> {
    pub fn new(refs: &'a RefSeqCollection) -> Self {
        Self { ref_id: 0, refs }
    }
}

impl<'a> Iterator for RefSeqIter<'a> {
    type Item = RefSeqSlice<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.ref_id < self.refs.n_refs() {
            let ref_id = self.ref_id;
            let slice = self.refs.get_refseq(ref_id);

            self.ref_id += 1;
            Some(slice)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {
    use crate::pf1::{cpp::DeserializeFromCpp, dense_index::DenseIndex, test_utils::*};

    use super::RefSeqCollection;

    #[test]
    fn test_tiny_from_fa() {
        let p = to_abs_path(TINY_REFS_FASTA);
        let refs = RefSeqCollection::from_fasta(p).unwrap();
        let r0 = "AGTGATGATAGTAGAGGTA";
        let r1 = "AGTGACTGATAGTAGCAGGTA";
        assert_eq!(refs.get_refseq(0).to_string(), r0);
        assert_eq!(refs.get_refseq(1).to_string(), r1);
    }

    #[test]
    fn test_refseq_kmer_iter_from_serialized_index() {
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(p).unwrap();

        let ref1 = pi.as_refseqs().get_refseq(1);

        assert_eq!(ref1.len(), 21);

        let k = pi.k();

        let kmers = vec![
            "agtga", "gtgac", "tgact", "gactg", "actga", "ctgat", "tgata", "gatag", "atagt",
            "tagta", "agtag", "gtagc", "tagca", "agcag", "gcagg", "caggt", "aggta",
        ];

        let kmers_ = Vec::from_iter(ref1.iter_kmers(k));
        let kmers_: Vec<String> = kmers_.iter().map(|km| format!("{}", km)).collect();
        assert_eq!(kmers, kmers_);
    }

    #[test]
    fn test_contig_iter() {
        // Check unitig tiling of small example
        // See: /test_data/tiny-multi-refs
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(p).unwrap();

        let refs = Vec::from_iter(pi.iter_refs());

        let ctgs0 = Vec::from_iter(pi.iter_unitigs_on_ref(&refs[0]));
        let ctgs1 = Vec::from_iter(pi.iter_unitigs_on_ref(&refs[1]));

        let lens0: Vec<usize> = ctgs0.iter().map(|x| x.unitig_len).collect();
        // let occs0: Vec<usize> = ctgs0.iter().map(|x| x.num_occs).collect();
        let ids0: Vec<usize> = ctgs0.iter().map(|x| x.unitig_id).collect();

        let lens1: Vec<usize> = ctgs1.iter().map(|x| x.unitig_len).collect();
        // let occs1: Vec<usize> = ctgs1.iter().map(|x| x.num_occs).collect();
        let ids1: Vec<usize> = ctgs1.iter().map(|x| x.unitig_id).collect();

        assert_eq!(ctgs0.len(), 5);
        assert_eq!(lens0, vec![5, 8, 9, 8, 5]);
        // assert_eq!(occs0, vec![2, 1, 2, 1, 2]);
        assert_eq!(ids0, vec![0, 1, 2, 3, 4]);

        assert_eq!(ctgs1.len(), 5);
        assert_eq!(lens1, vec![5, 9, 9, 9, 5]);
        // assert_eq!(occs1, vec![2, 1, 2, 1, 2]);
        assert_eq!(ids1, vec![0, 5, 2, 6, 4]);
    }
}
