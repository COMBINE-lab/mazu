use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use kmers::naive_impl::{CanonicalKmer, CanonicalKmerIterator, MatchType};

use crate::{
    kphf::{K2UPos, K2U},
    unitig_set::UnitigSet,
    MappedRefPos, ModIndex, ProjectedHits, U2Pos, UnitigOcc,
};

pub struct StreamingK2U<'a, T> {
    is_warm: bool,
    prev_k2upos: K2UPos,
    k2u: &'a T,
}

pub struct CachingU2Pos<'a, T>
where
    T: U2Pos,
{
    is_warm: bool,
    prev_uid: usize,
    prev_encoded_occs: T::EncodedOccs<'a>,

    #[allow(dead_code)]
    decode_is_warm: bool,
    #[allow(dead_code)]
    prev_decoded_uid: usize,
    #[allow(dead_code)]
    prev_decoded_occs: Vec<UnitigOcc>,
    u2pos: &'a T,
}

pub struct StreamingIndex<'a, K, U: U2Pos> {
    k2u: StreamingK2U<'a, K>,
    u2pos: CachingU2Pos<'a, U>,
}

impl<'a, T> StreamingK2U<'a, T>
where
    T: K2U + AsRef<UnitigSet>,
{
    pub fn new(k2u: &'a T) -> Self {
        let prev_k2upos = K2UPos {
            unitig_id: usize::MAX,
            unitig_len: usize::MAX,
            pos: usize::MAX,
            o: MatchType::NoMatch,
        };
        Self {
            is_warm: false,
            prev_k2upos,
            k2u,
        }
    }

    pub fn k(&self) -> usize {
        self.k2u.k()
    }

    pub fn k2u_streaming(&mut self, km: &CanonicalKmer) -> Option<K2UPos> {
        if self.is_warm {
            self.k2u_warm(km)
        } else {
            self.k2u_cold(km)
        }
    }

    fn k2u_warm(&mut self, km: &CanonicalKmer) -> Option<K2UPos> {
        let unitigs = self.k2u.as_ref();
        let k = unitigs.k();

        let next_pos = self.prev_k2upos.pos + 1;
        let last_km_pos = self.prev_k2upos.unitig_len - k;

        if next_pos > last_km_pos {
            // Not next km pos is not valid.
            return self.k2u_cold(km);
        }

        let prev_ui = self.prev_k2upos.unitig_id;
        let kw = unitigs.unitig_seq(prev_ui).get_kmer(next_pos, k);
        let mt = km.get_kmer_equivalency(&kw);

        match mt {
            MatchType::NoMatch => self.k2u_cold(km),
            _ => {
                self.prev_k2upos.pos = next_pos;
                self.prev_k2upos.o = mt;
                Some(self.prev_k2upos.clone())
            }
        }
    }

    fn k2u_cold(&mut self, km: &CanonicalKmer) -> Option<K2UPos> {
        self.prev_k2upos = self.k2u.k2u(km)?;
        self.is_warm = true;
        Some(self.prev_k2upos.clone())
    }
}

impl<'a, T> CachingU2Pos<'a, T>
where
    T: U2Pos,
    T::EncodedOccs<'a>: Clone,
{
    pub fn new(u2pos: &'a T) -> Self {
        Self {
            is_warm: false,
            prev_uid: usize::MAX,
            prev_encoded_occs: u2pos.empty_encoded_occs(),
            decode_is_warm: false,
            prev_decoded_uid: usize::MAX,
            prev_decoded_occs: Vec::new(),
            u2pos,
        }
    }
    pub fn encoded_unitig_occs(&mut self, uid: usize) -> T::EncodedOccs<'a> {
        if self.is_warm && uid == self.prev_uid {
            self.prev_encoded_occs.clone()
        } else {
            self.is_warm = true;
            self.prev_encoded_occs = self.u2pos.encoded_unitig_occs(uid);
            self.prev_encoded_occs.clone()
        }
    }

    pub fn decode_unitig_occs<MetaT>(
        &self,
        meta: &MetaT,
        occs: T::EncodedOccs<'a>,
    ) -> Vec<UnitigOcc> {
        // TODO cache decoded results.
        self.u2pos.decode_unitig_occs(meta, occs)
    }
}

impl<'a, K, U> StreamingIndex<'a, K, U>
where
    U: U2Pos,
    U::EncodedOccs<'a>: Clone,
    K: K2U + AsRef<UnitigSet>,
{
    pub fn k(&self) -> usize {
        self.k2u.k()
    }

    pub fn get_ref_pos(&mut self, km: &CanonicalKmer) -> Option<ProjectedHits<U::EncodedOccs<'a>>> {
        let k2upos = self.k2u.k2u_streaming(km)?;
        let uid = k2upos.unitig_id;
        let occs = self.u2pos.encoded_unitig_occs(uid);

        Some(ProjectedHits {
            kmer: km.clone(),
            k2upos,
            encoded_unitig_occs: occs,
        })
    }

    fn get_ref_pos_eager(&mut self, km: &CanonicalKmer) -> Option<Vec<MappedRefPos>> {
        let hits = self.get_ref_pos(km)?;
        Some(self.project_hits(hits))
    }

    pub fn project_hits(&self, hits: ProjectedHits<U::EncodedOccs<'a>>) -> Vec<MappedRefPos> {
        let u_occs = self
            .u2pos
            .decode_unitig_occs(&self.k2u, hits.encoded_unitig_occs);
        crate::index::project_onto_u_occs(self.k(), hits.k2upos, &u_occs)
    }

    // Validate canonical k-mers of u8 slice with given reference ID
    pub fn validate_ckmers(&mut self, ref_id: usize, seq: &[u8]) {
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
    pub fn validate_fasta<P: AsRef<Path>>(&mut self, fp: P) {
        let f = File::open(fp).unwrap();
        let f = BufReader::new(f);

        let mut lines = f.lines();
        let _k = self.k();

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

impl<'a, K, U> ModIndex<K, U>
where
    U: U2Pos,
    U::EncodedOccs<'a>: Clone,
    K: K2U + AsRef<UnitigSet>,
    Self: 'a,
{
    pub fn as_streaming(&'a self) -> StreamingIndex<K, U> {
        let k2u = StreamingK2U::new(self.as_k2u());
        let u2pos = CachingU2Pos::new(self.as_u2pos());
        StreamingIndex { k2u, u2pos }
    }
}

#[cfg(test)]
mod test {
    use crate::index::defaults::PiscemIndexDefault as PiscemIndex;
    use crate::test_utils::*;
    #[test]
    fn streaming_piscem_yeast() {
        let w = 15;
        let skew_param = 32;
        let prefix = YEAST_CF_PREFIX;
        let fasta = YEAST_CF_FA;

        let index = PiscemIndex::from_cf_prefix(prefix, w, skew_param).unwrap();

        // just some sanity checks
        assert!(index.as_sshash().n_minimizers() < index.n_kmers());
        assert_ne!(index.as_sshash().n_kmers_in_skew_index(), 0); // make sure the test case is not trivial

        let mut index = index.as_streaming();
        index.validate_fasta(fasta);
    }
}
