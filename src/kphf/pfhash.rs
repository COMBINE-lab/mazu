use kmers::naive_impl::{seq_vector::SeqVectorSlice, CanonicalKmer, MatchType};
use serde::{Deserialize, Serialize};
use simple_sds::{
    bit_vector::BitVector,
    int_vector::IntVector,
    ops::{Access, Pack},
    ops::{BitVec, Rank},
};

use rayon::prelude::*;

use super::{K2UPos, K2U, MPHF};
use crate::km_size_t;
use crate::unitig_set::UnitigSet;

#[allow(non_camel_case_types)]
type pf1_mphf_t = boomphf::Mphf<u64>;
pub type PFHashDefault = PFHash<pf1_mphf_t>;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct PFHash<MPHF> {
    unitigs: UnitigSet,
    mphf: MPHF,
    // pos: Vec<usize>,
    pos: IntVector,
}

impl<T> AsRef<UnitigSet> for PFHash<T> {
    fn as_ref(&self) -> &UnitigSet {
        &self.unitigs
    }
}

impl<T> PFHash<T> {
    pub fn from_parts(unitigs: UnitigSet, mphf: T, pos: IntVector) -> Self {
        Self { unitigs, mphf, pos }
    }
}

impl PFHash<pf1_mphf_t> {
    pub fn from_unitig_set(unitigs: UnitigSet) -> Self {
        let iter_kmers = unitigs.chunked_unitig_canonical_kmers();
        log::info!("Building BooMPHF");
        let mphf = pf1_mphf_t::from_chunked_iterator_parallel(
            1.7,
            &iter_kmers,
            None,
            unitigs.n_kmers(),
            rayon::current_num_threads(),
        );
        let mut pos = vec![usize::MAX; unitigs.n_kmers()];
        let ptr = crate::util::UnsafeSlice::new(&mut pos);

        log::info!("Inserting {} Kmer positions", unitigs.n_kmers());

        let offsets: Vec<usize> = (0..unitigs.n_unitigs())
            .map(|ui| unitigs.unitig_len(ui))
            .collect();
        let offsets = crate::util::prefix_sum(&offsets);

        (0..unitigs.n_unitigs()).into_par_iter().for_each(|ui| {
            for (km_pos, km) in unitigs.unitig_seq(ui).iter_kmers(unitigs.k()).enumerate() {
                let h = mphf.hash(&km.to_canonical().into_u64()) as usize;
                // pos[h] = km_pos + offsets[ui];
                unsafe {
                    ptr.write(h, km_pos + offsets[ui]);
                }
            }
        });

        let mut pos = IntVector::from(pos);
        pos.pack();
        Self { mphf, unitigs, pos }
    }
}

impl MPHF for pf1_mphf_t {
    fn try_hash_u64(&self, w: u64) -> Option<usize> {
        let h = self.try_hash(&w)?;
        Some(h as usize)
    }
}

impl<H: MPHF> K2U for PFHash<H> {
    fn k(&self) -> km_size_t {
        self.unitigs.k()
    }

    fn n_kmers(&self) -> usize {
        self.unitigs.n_kmers()
    }

    fn sum_unitigs_len(&self) -> usize {
        self.unitigs.total_len()
    }

    fn unitig_len(&self, id: usize) -> usize {
        self.unitigs.unitig_len(id)
    }

    fn n_unitigs(&self) -> usize {
        self.unitigs.n_unitigs()
    }

    fn unitig_seq(&self, id: usize) -> SeqVectorSlice {
        self.unitigs.unitig_seq(id)
    }

    fn k2u(&self, km: &CanonicalKmer) -> Option<K2UPos> {
        assert_eq!({ km.len() }, self.k());

        let word = km.get_canonical_word();
        let h = self.mphf.try_hash_u64(word)?;
        let km_pos = self.pos.get(h) as usize;

        let kw = self.unitigs.get_kmer_from_useq_pos(km_pos);

        let mt = km.get_kmer_equivalency(&kw);

        if let MatchType::NoMatch = mt {
            None
        } else {
            let unitig_id = self.unitigs.pos_to_id(km_pos);
            let unitig_len = self.unitig_len(unitig_id);
            let pos = km_pos - self.unitigs.unitig_start_pos(unitig_id);

            let ans = K2UPos {
                unitig_id,
                unitig_len,
                pos,
                o: mt,
            };
            Some(ans)
        }
    }
}

// Supports load-only compatibility with pufferfish 1 C++ implementation
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SampledPFHash<MPHF> {
    pub(crate) unitigs: UnitigSet,
    pub(crate) mphf: MPHF,

    pub(crate) sampled_pos: IntVector,   // sampled_positions
    pub(crate) sampled_vec: BitVector,   // is position sampled
    pub(crate) canonical_vec: BitVector, // orientation for non-sampled kmers
    pub(crate) direction_vec: BitVector, // which direction to walk
    pub(crate) ext_sizes: IntVector,     // extension sizes
    pub(crate) ext_bases: IntVector,     // the extension nucleotides
    pub(crate) sample_size: usize,
    pub(crate) extension_size: usize,
}

impl<H> AsRef<UnitigSet> for SampledPFHash<H> {
    fn as_ref(&self) -> &UnitigSet {
        &self.unitigs
    }
}

impl SampledPFHash<pf1_mphf_t> {
    pub fn from_unitig_set(_unitigs: UnitigSet) -> Self {
        todo!()
    }
}

impl<H: MPHF> K2U for SampledPFHash<H> {
    fn k(&self) -> km_size_t {
        self.unitigs.k()
    }

    fn n_kmers(&self) -> usize {
        self.unitigs.n_kmers()
    }

    fn sum_unitigs_len(&self) -> usize {
        self.unitigs.total_len()
    }

    fn unitig_len(&self, id: usize) -> usize {
        self.unitigs.unitig_len(id)
    }

    fn n_unitigs(&self) -> usize {
        self.unitigs.n_unitigs()
    }

    fn unitig_seq(&self, id: usize) -> SeqVectorSlice {
        self.unitigs.unitig_seq(id)
    }

    fn k2u(&self, kmer: &CanonicalKmer) -> Option<K2UPos> {
        let kw = kmer.get_canonical_word();

        let idx = self.mphf.try_hash_u64(kw)?;

        if self.sampled_vec.get(idx) {
            let pos_i = self.sampled_vec.rank(idx);
            let pos = self.sampled_pos.get(pos_i) as usize;
            self.k2u_w_pos(kmer, pos)
        } else {
            // get the pos.
            let mut signed_shift = 0i64;

            let current_rank = self.sampled_vec.rank(idx); // h(x) is not sampled
            let extension_pos = idx - current_rank; // get position of extension, minus curr_rank since curr_rank sampled kmers have no extensions
            let extension_word = self.ext_bases.get(extension_pos); // get extension bases.

            // shift to find the nearest sampled `km`
            let sampled_km = {
                let mut km = kmer.clone();

                // if fw is canonical and canonical vec is not canonical, swap
                // if fw is not canonical and canonical vec is canonical, swap. So XOR...
                if !self.canonical_vec.get(extension_pos) ^ !km.is_fw_canonical() {
                    km.swap();
                }

                let shift_fw = self.direction_vec.get(extension_pos);
                let llimit =
                    self.extension_size - (self.ext_sizes.get(extension_pos) as usize + 1usize);

                if shift_fw {
                    for i in (llimit + 1..=self.extension_size).rev() {
                        let ssize = 2 * (i - 1);
                        let curr_code = (extension_word & (0x3 << ssize)) >> ssize;
                        km.append_base(curr_code);
                        signed_shift -= 1;
                    }
                } else {
                    for i in (llimit + 1..=self.extension_size).rev() {
                        let ssize = 2 * (i - 1);
                        let curr_code = (extension_word & (0x3 << ssize)) >> ssize;
                        km.prepend_base(curr_code);
                        signed_shift += 1;
                    }
                }
                km
            };

            // compute pos of sampled km
            let kw = sampled_km.get_canonical_word();
            let idx = self.mphf.try_hash_u64(kw)?;
            if !self.sampled_vec.get(idx) {
                return None;
            }

            let current_rank = self.sampled_vec.rank(idx);
            let sample_pos = self.sampled_pos.get(current_rank);

            // pos of kmer
            let pos = ((sample_pos as i64) + signed_shift) as usize;

            // make sure new pos does not cross unitig boundary
            if self.unitigs.is_valid_useq_pos(pos) {
                self.k2u_w_pos(kmer, pos)
            } else {
                None
            }
        }
    }
}

impl<H: MPHF> SampledPFHash<H> {
    fn k2u_w_pos(&self, km: &CanonicalKmer, pos: usize) -> Option<K2UPos> {
        // return k2u given candidate pos on useq
        let kw = self.unitigs.get_kmer_from_useq_pos(pos);

        let mt = km.get_kmer_equivalency(&kw);

        if let MatchType::NoMatch = mt {
            None
        } else {
            let unitig_id = self.unitigs.pos_to_id(pos);
            let unitig_len = self.unitig_len(unitig_id);
            let pos = pos - self.unitigs.unitig_start_pos(unitig_id);

            let ans = K2UPos {
                unitig_id,
                unitig_len,
                pos,
                o: mt,
            };
            Some(ans)
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::test_utils::*;

    #[test]
    fn tiny() {
        // tiny.seq:
        // 3	CACACACCAC
        // 2	CCTCAATACG
        let unitigs = load_unitigs(TINY_CF_PREFIX);
        let pfhash = PFHash::from_unitig_set(unitigs);
        pfhash.validate_self()
    }

    #[test]
    #[should_panic]
    fn tiny_kmer_too_small() {
        let unitigs = load_unitigs(TINY_CF_PREFIX);
        let pfhash = PFHash::from_unitig_set(unitigs);
        let km = CanonicalKmer::from_u64(0, (pfhash.k() + 1) as u8);
        pfhash.k2u(&km);
    }
}
