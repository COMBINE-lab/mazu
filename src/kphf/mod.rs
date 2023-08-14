use kmers::naive_impl::{seq_vector::SeqVectorSlice, CanonicalKmer, MatchType};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::hash::BuildHasher;
use wyhash::WyHash;

pub mod pfhash;
pub mod sshash;

pub use pfhash::{PFHash, SampledPFHash};
pub use sshash::SSHash;

#[derive(Debug, Clone, PartialEq)]
pub struct K2UPos {
    pub unitig_id: usize,
    pub unitig_len: usize,
    pub pos: usize,
    pub o: MatchType,
}

impl K2UPos {
    /// Reverse the match type (like `std::cmp::Ordering::reverse()`)
    pub fn reverse_match_type(&mut self) {
        self.o = match self.o {
            MatchType::IdentityMatch => MatchType::TwinMatch,
            MatchType::TwinMatch => MatchType::IdentityMatch,
            MatchType::NoMatch => MatchType::NoMatch,
        };
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct WyHashState(u64);

impl WyHashState {
    pub fn with_seed(seed: u64) -> Self {
        Self(seed)
    }
}

impl Default for WyHashState {
    fn default() -> Self {
        Self::with_seed(0)
    }
}

impl BuildHasher for WyHashState {
    type Hasher = WyHash;
    fn build_hasher(&self) -> Self::Hasher {
        WyHash::with_seed(self.0)
    }
}

pub trait MPHF {
    fn try_hash_u64(&self, _w: u64) -> Option<usize>;
}

pub trait K2U {
    fn k(&self) -> usize;
    // fn k2u_fw(&self, km: &Kmer) -> Option<K2UPos>;
    fn unitig_len(&self, id: usize) -> usize;
    fn n_unitigs(&self) -> usize;
    fn unitig_seq(&self, id: usize) -> SeqVectorSlice;
    fn n_kmers(&self) -> usize;
    fn sum_unitigs_len(&self) -> usize;
    fn k2u(&self, km: &CanonicalKmer) -> Option<K2UPos>;

    // Serial self validation for sane debugging
    fn validate_self(&self)
    where
        Self: Sync,
    {
        (0..self.n_unitigs()).for_each(|ui| {
            let useq = self.unitig_seq(ui);
            for (pos, km) in useq.iter_kmers(self.k()).enumerate() {
                let unitig_len = self.unitig_len(ui);

                let mut km = CanonicalKmer::from(km);
                let ans = K2UPos {
                    pos,
                    unitig_len,
                    unitig_id: ui,
                    o: MatchType::IdentityMatch,
                };
                let res = self.k2u(&km);
                if res.is_none() {
                    panic!("{} not found at pos {} on unitig {}", km, pos, ui)
                }
                let res = res.unwrap();
                assert_eq!(ans, res);

                km.swap();
                let ans = K2UPos {
                    pos,
                    unitig_len,
                    unitig_id: ui,
                    o: MatchType::TwinMatch,
                };
                let res = self.k2u(&km).unwrap();
                assert_eq!(ans, res);
            }
        })
    }

    fn validate_self_parallel(&self)
    where
        Self: Sync,
    {
        (0..self.n_unitigs()).into_par_iter().for_each(|ui| {
            let useq = self.unitig_seq(ui);
            for (pos, km) in useq.iter_kmers(self.k()).enumerate() {
                let unitig_len = self.unitig_len(ui);

                let mut km = CanonicalKmer::from(km);
                let ans = K2UPos {
                    pos,
                    unitig_len,
                    unitig_id: ui,
                    o: MatchType::IdentityMatch,
                };
                let res = self.k2u(&km);
                if res.is_none() {
                    panic!("{} not found at pos {} on unitig {}", km, pos, ui)
                }
                let res = res.unwrap();
                assert_eq!(ans, res);

                km.swap();
                let ans = K2UPos {
                    pos,
                    unitig_len,
                    unitig_id: ui,
                    o: MatchType::TwinMatch,
                };
                let res = self.k2u(&km).unwrap();
                assert_eq!(ans, res);
            }
        });
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::test_utils::*;

    #[test]
    fn yeast_sshash() {
        let unitigs = load_unitigs(YEAST_CF_PREFIX);
        let w = 15;
        let bh = WyHashState::default();
        let sshash = SSHash::from_unitig_set_no_skew_index(unitigs, w, bh).unwrap();
        sshash.validate_self()
    }

    #[test]
    fn yeast_pfhash() {
        let unitigs = load_unitigs(YEAST_CF_PREFIX);
        let pfhash = PFHash::from_unitig_set(unitigs);
        pfhash.validate_self()
    }
}
