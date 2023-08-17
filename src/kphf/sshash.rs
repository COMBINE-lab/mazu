use kmers::naive_impl::{seq_vector::SeqVectorSlice, CanonicalKmer, Kmer, MatchType};
use log;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use simple_sds::{
    int_vector::IntVector,
    ops::{Access, Pack, Vector},
};

use std::hash::BuildHasher;

use super::{K2UPos, WyHashState, K2U};
use crate::{elias_fano::EFVector, unitig_set::UnitigSet, Result};

#[allow(non_camel_case_types)]
pub(crate) type mphf_t = boomphf::Mphf<u64>;
pub type SSHashDefault = SSHash<mphf_t, WyHashState>;
pub struct SSHashBuilder<MPHF, T: BuildHasher> {
    w: usize,
    unitigs: UnitigSet,
    mphf: MPHF,
    build_hasher: T,
    occs_prefix_sum: Vec<usize>,
    pos: Vec<usize>,

    skew_param: usize,
    skew_index: Option<SkewIndex<MPHF>>,
}

// An SSHash implementation with some differences:
// - Optimized to skip linear scan of super-k-mer positions.
// - Uses offset of minimizer in queried k-mer to immediately retrieve candidate fw and rc k-mer positions
// - achieved by defining mini(g*) = mini(min(g, g')) instead of min(mini(g), mini(g'))
// - A flat, single level skew index k-mers directly to positins, instead of k-mers to minimizer positions.
//   - FIXME: could update to point to minimizer positions instead of k-mer positions.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SSHash<MPHF, T: BuildHasher> {
    build_hasher: T, // buildhasher for minimizer computation

    // k: usize,
    w: usize,
    unitigs: UnitigSet,
    mphf: MPHF,
    occs_prefix_sum: EFVector, // the "Sizes" of buckets (as in paper)
    pos: IntVector,            // the "Offsets" of minimizers (as in paper)

    skew_param: usize,
    skew_index: Option<SkewIndex<MPHF>>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SkewIndex<MPHF> {
    mphf: MPHF,
    pos: IntVector, // k-mer positions
}

impl SkewIndex<mphf_t> {
    #[inline]
    fn lookup(&self, km: &CanonicalKmer) -> Option<usize> {
        let word = km.get_canonical_word();
        let h = self.mphf.try_hash(&word)? as usize;
        let km_pos = self.pos.get(h) as usize;
        Some(km_pos)
    }

    fn n_kmers(&self) -> usize {
        self.pos.len()
    }

    fn num_bits(&self) -> usize {
        log::warn!("External MPHF implementation does not support 'num_bits()', approximating with bincode::serialized_size() instead...");
        let mphf_bytes = bincode::serialized_size(&self.mphf).unwrap() as usize;
        let mphf_bits = mphf_bytes * 8;
        self.pos.num_bits() + mphf_bits
    }
}

impl<BH: BuildHasher + Clone> SSHashBuilder<mphf_t, BH> {
    pub fn from_unitig_set_no_skew_index(unitigs: UnitigSet, w: usize, build_hasher: BH) -> Self {
        let skew_param = usize::MAX;
        Self::from_unitig_set(unitigs, w, skew_param, build_hasher)
    }

    pub fn from_unitig_set(
        unitigs: UnitigSet,
        w: usize,
        skew_param: usize,
        build_hasher: BH,
    ) -> Self {
        assert!(w <= unitigs.k());

        // 1. collect minimizers
        log::info!("Collecting all canonical super-kmers (and corresponding minimizer positions)");

        let mut minimizers = Vec::new();
        let mut useq_offset = 0;

        for ui in 0..unitigs.n_unitigs() {
            let u = unitigs.unitig_seq(ui);

            {
                let fw_mm_iter = u.iter_canonical_minimizers(unitigs.k(), w, build_hasher.clone())
                    .filter(|(mm, is_fw_canonical)| *is_fw_canonical)
                    .map(|(mm, _)| mm);

                let mut prev_mm = None;
                for mut mm in fw_mm_iter {
                    mm.pos += useq_offset;
                    let curr_mm = Some(mm.clone());
                    if curr_mm != prev_mm {
                        minimizers.push(mm.clone());
                    }
                    prev_mm = curr_mm;
                }
            }

            {
                let rc_mm_iter = u.iter_canonical_minimizers(unitigs.k(), w, build_hasher.clone())
                    .filter(|(mm, is_fw_canonical)| !*is_fw_canonical)
                    .map(|(mm, _)| mm);

                let mut prev_mm = None;
                for mut mm in rc_mm_iter {
                    mm.pos += useq_offset;
                    let curr_mm = Some(mm.clone());
                    if curr_mm != prev_mm {
                        minimizers.push(mm.clone());
                    }
                    prev_mm = curr_mm;
                }
            }

            // for mut skmer in u.iter_canonical_super_kmers(unitigs.k(), w, build_hasher.clone()) {
            //     skmer.inc_pos(useq_offset);

            //     super_kmers.push(skmer);
            // }
            useq_offset += u.len();
        }

        // 2. Sort and deduplicate
        log::info!("Found {} unique minimizer occs", minimizers.len());
        log::info!("Deduplicating minimizer occurrences by sorting... ");

        log::info!("\t* sorting unique minimizer occurrences...");
        minimizers.par_sort_by_key(|sk| sk.as_u64()); // sort occurrences by their sequence

        log::info!("\t* deduplicating to get minimizer set...");
        let mut mm_occs = Vec::new(); // mm_occs[i] is the number of times minimizer_i occurs
        let mut mm_set = Vec::new(); // mm_set[i] contains minimizer_i
        {
            let mut curr_mm = minimizers[0].as_u64();
            let mut curr_occs = 0_usize;

            for mm in &minimizers {
                // pos.push(mm.pos);
                if curr_mm != mm.as_u64() {
                    mm_occs.push(curr_occs);
                    mm_set.push(curr_mm);

                    curr_mm = mm.as_u64();
                    curr_occs = 0;
                }
                curr_occs += 1;
            }
            mm_occs.push(curr_occs);
            mm_set.push(curr_mm);
        }

        log::info!("Constructing BooMPHF... ");

        // 3) Build MPHF over minmizers
        let mphf = mphf_t::new_parallel(1.7, &mm_set, None);

        // 4) Build prefix sum
        log::info!("Building prefix sum of minimizer counts... ");
        let occs_prefix_sum;
        {
            let mut n_occs = vec![usize::MAX; mm_set.len()];
            for (mm, &n) in mm_set.iter().zip(&mm_occs) {
                let h = mphf.hash(mm) as usize;
                n_occs[h] = n;
            }
            occs_prefix_sum = crate::util::prefix_sum(&n_occs);
        }

        // 5) Insert Minimizer positions according to mphf values
        log::info!("Inserting minimizer positions... ");
        log::info!("\t* {} unique minimizers", mm_set.len());
        log::info!("\t* {} total minimizer occurrences", minimizers.len());

        let mut pos = vec![usize::MAX; minimizers.len()];
        let ptr: crate::util::UnsafeSlice<'_, usize> = crate::util::UnsafeSlice::new(&mut pos);

        // prefix sum over number of occurrences of mmer_i in mm_occs[i] yields range
        // in vec of super-k-mers containing corresponding super-k-mers
        let ranges = crate::util::prefix_sum(&mm_occs);
        let _ = &mm_set.par_iter().enumerate().for_each(|(i, mm)| {
            // 1) hash the value of mmer in minimizer set
            let h = mphf.hash(mm) as usize;

            // 2) get index of interval pointing at corresponding super k-mers
            let s = ranges[i];
            let e = ranges[i + 1];

            // 3) get index of interval corresponding to hash value
            let sh = occs_prefix_sum[h];

            let positions = minimizers[s..e].iter().map(|mm| mm.pos);
            for (i, p) in positions.enumerate() {
                unsafe {
                    ptr.write(sh + i, p);
                }
            }
        });

        // 6) build skew index
        let skew_index = if skew_param == usize::MAX {
            None
        } else {
            log::info!("Building skew index");
            log::info!("\t* collecting skew super-k-mers...");
            // let mut km_set = Vec::new();
            // let mut km_positions = Vec::new();
            let mut skew_tuples = Vec::new();
            let skew_min = skew_param;
            for (i, &occs) in mm_occs.iter().enumerate() {
                if occs <= skew_min {
                    continue;
                }
                // otherwise stick it in the skew index!

                // Get indexes of interval pointing at corresponding minimizers
                let s = ranges[i];
                let e = ranges[i + 1];
                for mm in &minimizers[s..e] {
                    let start_pos = {
                        if mm.pos < (unitigs.k() - w) {
                            0
                        } else {
                            mm.pos - (unitigs.k() - w)
                        }
                    };
                    let n_kmers = unitigs.k() - w + 1;

                    // let start_pos = sk.start_pos();
                    // let n_kmers = sk.n_kmers();
                    // collect kmers overlapping minimizer positionss
                    for offset in 0..n_kmers {

                        let pos = start_pos + offset;
                        if unitigs.is_valid_useq_pos(pos) {
                            let km = unitigs.get_kmer_from_useq_pos(pos);
                            let word = km.to_canonical_word();
                            skew_tuples.push((word, pos));
                        }
                    }
                }
            }

            // deduplicate skew tuples
            log::info!("\t* extracted {} kmers overlapping skew minimizers", skew_tuples.len());
            log::info!("\t* deduplicating {} candidate skew kmers", skew_tuples.len());
            skew_tuples.par_sort_by_key(|tup| tup.0);
            skew_tuples.dedup_by_key(|tup| tup.0);
            let km_set: Vec<u64> = skew_tuples.iter().map(|tup| tup.0).collect();

            log::info!("\t* indexing {} kmers in skew index...", km_set.len());
            log::info!("\t* building skew mphf...");

            let skew_mphf = mphf_t::new_parallel(1.7, &km_set, None);
            let mut skew_pos = vec![0; km_set.len()];

            log::info!("\t* inserting skew k-mer positions");
            for (km, pos) in skew_tuples {
                let h = skew_mphf.try_hash(&km).unwrap() as usize;
                skew_pos[h] = pos;
            }

            log::info!("\t* packing...");
            let mut skew_pos = IntVector::from(skew_pos);
            skew_pos.pack();
            Some(SkewIndex {
                mphf: skew_mphf,
                pos: skew_pos,
            })
        };

        Self {
            w,
            unitigs,
            mphf,
            build_hasher,
            occs_prefix_sum,
            pos,
            skew_index,
            skew_param,
        }
    }

    pub fn finish(self) -> Result<SSHash<mphf_t, BH>> {
        log::info!("Elias Fano encoding int vector of positions");
        let occs_prefix_sum = EFVector::from_usize_slice(&self.occs_prefix_sum)?;

        let mut pos = IntVector::from(self.pos);
        pos.pack();

        log::info!("Done.");
        Ok(SSHash {
            // k: self.k,
            w: self.w,
            unitigs: self.unitigs,
            mphf: self.mphf,
            build_hasher: self.build_hasher,
            occs_prefix_sum,
            pos,
            skew_index: self.skew_index,
            skew_param: self.skew_param,
        })
    }
}

impl<T: BuildHasher + Clone> SSHash<mphf_t, T> {
    pub fn n_minimizers(&self) -> usize {
        self.occs_prefix_sum.len()
    }

    pub fn n_kmers_in_skew_index(&self) -> usize {
        self.skew_index.as_ref().unwrap().n_kmers()
    }

    pub fn num_bits(&self) -> usize {
        // Note: https://docs.rs/bincode/1.3.3/src/bincode/lib.rs.html#193-201,
        // we don't have to worry about bincode changing integer widths under the hood
        log::warn!("External MPHF implementation does not support 'num_bits()', approximating with bincode::serialized_size() instead...");

        let mphf_bytes = bincode::serialized_size(&self.mphf).unwrap() as usize;
        let skew_index = match &self.skew_index {
            None => 0,
            Some(skew_index) => skew_index.num_bits(),
        };

        std::mem::size_of::<usize>() * 8
            + self.unitigs.num_bits()
            + self.occs_prefix_sum.num_bits()
            + self.pos.num_bits()
            + mphf_bytes * 8
            + skew_index
    }

    pub fn print_stats(&self) {
        log::info!("*** STATISTICS");
        log::info!("kmers: {}", self.n_kmers());
        log::info!("unitigs: {}", self.n_unitigs());
        log::info!("Total size:      {} bytes", self.num_bits() / 8);
        log::info!(
            "occs_prefix_sum: {} bytes",
            self.occs_prefix_sum.num_bits() / 8
        );

        dbg!(self.occs_prefix_sum.len());
        dbg!(self.occs_prefix_sum.low_bit_width());
        dbg!(self.occs_prefix_sum.num_high_bits());

        dbg!(self.pos.width());
        dbg!(self.pos.len());
        dbg!(self.n_minimizers());
        dbg!(self.n_kmers());

        log::info!("pos:             {} bytes", self.pos.num_bits() / 8);
        log::info!("unitig set:      {} bytes", self.unitigs.num_bits() / 8);

        log::info!(
            "occs_prefix_sum: {} bpk",
            self.occs_prefix_sum.num_bits() as f64 / self.n_kmers() as f64
        );
        log::info!(
            "pos:             {} bpk",
            self.pos.num_bits() as f64 / self.n_kmers() as f64
        );
        log::info!(
            "unitig set:      {} bpk",
            self.unitigs.num_bits() as f64 / self.n_kmers() as f64
        );
        log::info!(
            "bits / kmer:     {}",
            (self.num_bits() as f64) / (self.n_kmers() as f64)
        );

        log::info!("*** UnitigSet stats");
        self.unitigs.print_stats();
    }

    pub fn from_unitig_set_no_skew_index(
        unitigs: UnitigSet,
        w: usize,
        build_hasher: T,
    ) -> Result<Self> {
        SSHashBuilder::from_unitig_set_no_skew_index(unitigs, w, build_hasher).finish()
    }

    pub fn from_unitig_set(
        unitigs: UnitigSet,
        w: usize,
        skew_param: usize,
        build_hasher: T,
    ) -> Result<Self> {
        SSHashBuilder::from_unitig_set(unitigs, w, skew_param, build_hasher).finish()
    }

    // Do K-mer to unitig lookup in skew index
    fn k2u_skew_index(&self, km: &CanonicalKmer) -> Option<K2UPos> {
        let pos = self.skew_index.as_ref()?.lookup(km)?;
        let kw = self.unitigs.get_kmer_from_useq_pos(pos);
        let mt = km.get_kmer_equivalency(&kw);
        match mt {
            MatchType::NoMatch => None,
            _ => {
                let unitig_id = self.unitigs.pos_to_id(pos);
                let unitig_len = self.unitig_len(unitig_id);
                let pos = pos - self.unitigs.unitig_start_pos(unitig_id);
                Some(K2UPos {
                    unitig_id,
                    unitig_len,
                    pos,
                    o: mt,
                })
            }
        }
    }
}

impl<H, BH: BuildHasher> AsRef<UnitigSet> for SSHash<H, BH> {
    fn as_ref(&self) -> &UnitigSet {
        &self.unitigs
    }
}

impl<T: BuildHasher + Clone> K2U for SSHash<mphf_t, T> {
    fn k(&self) -> usize {
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
        let km_fw = km.get_fw_mer();
        assert_eq!(km.get_fw_mer().k as usize, self.k());

        // get canonical minimizer from forward km.
        let (mmer, offset) = km_fw.canonical_minimizer(self.w, &self.build_hasher);

        let h = self.mphf.try_hash(&mmer.into_u64());
        let h = h? as usize;

        // get candidate minimizer occ positions
        let pos_start = self.occs_prefix_sum.get(h);
        let pos_end = self.occs_prefix_sum.get(h + 1);

        // check wrt projected positions w.r.t in fw and rc offsets
        let n_occs = (pos_end - pos_start) as usize;
        if n_occs > self.skew_param {
            // in skew index;
            return self.k2u_skew_index(km);
        }

        let last_km_start_pos = self.unitigs.total_len() - self.k();
        let rc_offset = self.k() - offset - self.w;
        for pi in pos_start..pos_end {
            let mm_pos: usize = self.pos.get(pi as usize) as usize;

            // Check candidate kmer from useq w.r.t offsets of fw_mer of queried canonical kmer
            if (mm_pos >= offset) && (mm_pos - offset) <= last_km_start_pos {
                let km_pos = mm_pos - offset;
                let kw = self.unitigs.get_kmer_u64_from_useq_pos(km_pos);

                // Check queried canonical kmer equivalancy with fw km word on useq.
                let mt = km.get_word_equivalency(kw);
                if !matches!(mt, MatchType::NoMatch) {
                    let unitig_id = self.unitigs.pos_to_id(km_pos);
                    let unitig_len = self.unitig_len(unitig_id);
                    let pos = km_pos - self.unitigs.unitig_start_pos(unitig_id);

                    // Check that the position is valid, there is a degenerate case where
                    // an invalid kmer encoded in useq but crossing a unitig-unitig boundary can
                    //  (1) contain a minimizer, and (2) be valid kmer in a valid unitig
                    // So check that the pos of last kmer char is at most pos of last char on unitig
                    let end_pos = km_pos + self.k();
                    if end_pos <= self.unitigs.unitig_end_pos(unitig_id) {
                        let ans = K2UPos {
                            unitig_id,
                            unitig_len,
                            pos,
                            o: mt,
                        };
                        return Some(ans);
                    }
                }
            }

            // Check candidate kmer from useq w.r.t offsets of rc_mer of queried canonical kmer
            if (mm_pos >= rc_offset) && (mm_pos - rc_offset) <= last_km_start_pos {
                let km_pos = mm_pos - rc_offset;
                let kw = self.unitigs.get_kmer_u64_from_useq_pos(km_pos);

                // Check queried canonical kmer equivalancy with fw km word on useq.
                let mt = km.get_word_equivalency(kw);
                if !matches!(mt, MatchType::NoMatch) {
                    let unitig_id = self.unitigs.pos_to_id(km_pos);
                    let unitig_len = self.unitig_len(unitig_id);
                    let pos = km_pos - self.unitigs.unitig_start_pos(unitig_id);

                    // Same unitig-unitig boundary check.
                    let end_pos = km_pos + self.k();

                    if end_pos <= self.unitigs.unitig_end_pos(unitig_id) {
                        let ans = K2UPos {
                            unitig_id,
                            unitig_len,
                            pos,
                            o: mt,
                        };
                        return Some(ans);
                    }
                }
            }
        }
        None
    }
}

impl<T: BuildHasher + Clone> SSHash<mphf_t, T> {
    // Get queried k-mer  offset into unitig sequence in the fw orientation only.
    // Notes:
    // - does not check if reverse k-mer occurs on unitig sequence.
    // - this fn was kept from an older version that did not support canonical k-ner queries
    // - is simpler and kept for unit tests.
    pub fn k2u_fw(&self, km: &Kmer) -> Option<K2UPos> {
        assert_eq!(km.k as usize, self.k());

        // get canonical minimizer from queried Kmer.
        let (mmer, offset) = km.canonical_minimizer(self.w, &self.build_hasher);

        let h = self.mphf.try_hash(&mmer.into_u64());
        let h = h? as usize;

        // get candidate minimizer occ positions
        let pos_start = self.occs_prefix_sum.get(h);
        let pos_end = self.occs_prefix_sum.get(h + 1);

        // we only need to check fw orientation of queried k-mer for equivalency
        for pi in pos_start..pos_end {
            // Note:
            // - perhaps boundchecking for kmer retrieval should be done by UnitigSet so
            //   internal representation is not exposed. But this is fine for now...
            let mm_pos = self.pos.get(pi as usize) as usize;

            // 1) check that minimizer postion doesn't push km_pos to -ve
            if mm_pos < offset {
                continue;
            }

            // 2) check that kmer position doesn't fall of last unitig position.
            let km_pos = mm_pos - offset;
            let last_start_pos = self.unitigs.total_len() - self.k();
            if km_pos > last_start_pos {
                continue;
            }

            // Notes:
            // - maybe wrap with unsafe? because it can cross a unitig-unitig boundary
            // - maybe change get_kmer func to return option to avoid 2) check above
            //   and exposing internal representation.
            let kw = self.unitigs.get_kmer_u64_from_useq_pos(km_pos);

            if kw == km.into_u64() {
                let unitig_id = self.unitigs.pos_to_id(km_pos);
                let unitig_len = self.unitig_len(unitig_id);
                let pos = km_pos - self.unitigs.unitig_start_pos(unitig_id);

                // Check that the position is valid, there is a degenerate case where
                // an invalid kmer encoded in useq but crossing a unitig-unitig boundary can
                //  (1) contain a minimizer, and (2) be valid kmer in a valid unitig
                let end_pos = km_pos + self.k();
                if end_pos > self.unitigs.unitig_end_pos(unitig_id) {
                    continue;
                }

                let ans = K2UPos {
                    unitig_id,
                    unitig_len,
                    pos,
                    o: MatchType::IdentityMatch,
                };
                return Some(ans);
            }
        }
        None
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::kphf::WyHashState;
    use crate::test_utils::*;

    #[test]
    fn tiny_w3() {
        let w = 3;
        test_tiny_w_params(w, 0)
    }

    #[test]
    fn tiny_w5() {
        let w = 5;
        test_tiny_w_params(w, 0)
    }

    fn test_tiny_w_params(w: usize, _seed: u64) {
        // tiny.seg (k=7)
        // 3	CACACACCAC
        // 2	CCTCAATACG
        let unitigs = load_unitigs(TINY_CF_PREFIX);
        let bh = WyHashState::default();
        let sshash = SSHash::from_unitig_set_no_skew_index(unitigs, w, bh).unwrap();

        // 1) fw
        let mut km = CanonicalKmer::from("CACACAC");
        let pos: K2UPos = sshash.k2u_fw(&km.get_fw_mer()).unwrap();
        let pos_fw: K2UPos = sshash.k2u(&km).unwrap();

        assert_eq!(pos, pos_fw);

        let mut ans = K2UPos {
            unitig_id: 0,
            unitig_len: 10,
            pos: 0,
            o: MatchType::IdentityMatch,
        };

        assert_eq!(pos, ans);
        // 1) rc
        km.swap();
        let pos_rc: K2UPos = sshash.k2u(&km).unwrap();
        ans.reverse_match_type();
        assert_eq!(pos_rc, ans);

        // 2) fw
        let mut km = CanonicalKmer::from("ACACACC");
        let pos: K2UPos = sshash.k2u_fw(&km.get_fw_mer()).unwrap();
        let pos_fw: K2UPos = sshash.k2u(&km).unwrap();

        assert_eq!(pos, pos_fw);
        let mut ans = K2UPos {
            unitig_id: 0,
            unitig_len: 10,
            pos: 1,
            o: MatchType::IdentityMatch,
        };

        assert_eq!(pos, ans);
        km.swap();
        let pos_rc: K2UPos = sshash.k2u(&km).unwrap();
        ans.reverse_match_type();
        assert_eq!(pos_rc, ans);

        // 3)
        let mut km = CanonicalKmer::from("ACACCAC");
        let pos = sshash.k2u_fw(&km.get_fw_mer()).unwrap();
        let pos_fw: K2UPos = sshash.k2u(&km).unwrap();
        assert_eq!(pos, pos_fw);

        let mut ans = K2UPos {
            unitig_id: 0,
            unitig_len: 10,
            pos: 3,
            o: MatchType::IdentityMatch,
        };

        assert_eq!(pos, ans);

        km.swap(); // ACACCAC (gtggtgt) [tgtggtg]
        let pos_rc: K2UPos = sshash.k2u(&km).unwrap();
        ans.reverse_match_type();
        assert_eq!(pos_rc, ans);

        // 4)
        let mut km = CanonicalKmer::from("CCTCAAT");
        let pos = sshash.k2u_fw(&km.get_fw_mer()).unwrap();
        let pos_fw: K2UPos = sshash.k2u(&km).unwrap();
        assert_eq!(pos, pos_fw);

        let mut ans = K2UPos {
            unitig_id: 1,
            unitig_len: 10,
            pos: 0,
            o: MatchType::IdentityMatch,
        };
        assert_eq!(pos, ans);

        km.swap();
        let pos_rc: K2UPos = sshash.k2u(&km).unwrap();
        ans.reverse_match_type();
        assert_eq!(pos_rc, ans);

        // 5)
        // tiny.seg (k=7)
        // 3	CACACACCAC
        // 2	CCTCAATACG
        let mut km = CanonicalKmer::from("CAATACG");
        let pos = sshash.k2u_fw(&km.get_fw_mer()).unwrap();
        let pos_fw: K2UPos = sshash.k2u(&km).unwrap();
        assert_eq!(pos, pos_fw);

        let mut ans = K2UPos {
            unitig_id: 1,
            unitig_len: 10,
            pos: 3,
            o: MatchType::IdentityMatch,
        };

        assert_eq!(pos_fw, ans);
        km.swap();
        let pos_rc: K2UPos = sshash.k2u(&km).unwrap();
        ans.reverse_match_type();
        assert_eq!(pos_rc, ans);

        // Miss
        let km = Kmer::from("AAAAAAA");
        let pos = sshash.k2u_fw(&km);
        assert_eq!(pos, None);

        let mut km = CanonicalKmer::from(km);
        assert_eq!(sshash.k2u(&km), None);
        km.swap();
        assert_eq!(sshash.k2u(&km), None);
    }

    #[test]
    fn tiny_vary_seeds_and_windows() {
        let unitigs = load_unitigs(TINY_CF_PREFIX);

        for seed in 0..100 {
            for w in 1..=unitigs.k() {
                test_tiny_w_params(w, seed)
            }
        }
    }

    #[test]
    fn validate_self() {
        let unitigs = load_unitigs(TINY_CF_PREFIX);
        let w = 3;
        let bh = WyHashState::default();

        let sshash = SSHash::from_unitig_set_no_skew_index(unitigs.clone(), w, bh.clone()).unwrap();
        sshash.validate_self();

        let sshash = SSHash::from_unitig_set(unitigs, 5, 0, bh).unwrap();
        sshash.validate_self();
    }

    // TODO: write a unit test for canonical case where a minimizer occurrence induces a valid k-mer on unitig-unitig boundary
    // on concatenated unitig sequence.
    #[test]
    fn unitigs_share_mmer() {
        let seqs = vec![
            "ACAACTTACCCTCCATTACCCTACCTCCCCA".to_string(),
            "CAACTTACCCTCCATTACCCTACCTCCCCAC".to_string(),
        ];
        let w = 15;
        let bh = WyHashState::default();

        let unitigs = UnitigSet::from_seqs(&seqs, 31).unwrap();
        let sshash = SSHash::from_unitig_set_no_skew_index(unitigs, w, bh).unwrap();

        let k1 = Kmer::from(seqs[1].clone());
        let ans = K2UPos {
            pos: 0,
            unitig_len: seqs[1].len(),
            unitig_id: 1,
            o: MatchType::IdentityMatch,
        };
        assert_eq!(sshash.k2u_fw(&k1.to_reverse_complement()), None);
        assert_eq!(sshash.k2u_fw(&k1).unwrap(), ans);

        let k0 = Kmer::from(seqs[0].clone());
        let ans = K2UPos {
            pos: 0,
            unitig_len: seqs[0].len(),
            unitig_id: 0,
            o: MatchType::IdentityMatch,
        };
        assert_eq!(sshash.k2u_fw(&k0).unwrap(), ans);

        // todo!("This test fails, check trace, use lexhasher, and implement validate_self_serial for serial validation");
        sshash.validate_self()
    }

    #[test]
    #[should_panic]
    fn tiny_kmer_too_small() {
        let unitigs = load_unitigs(TINY_CF_PREFIX);
        let w = 3;
        let bh = WyHashState::default();

        let sshash = SSHash::from_unitig_set_no_skew_index(unitigs, w, bh).unwrap();
        sshash.validate_self();

        let km = Kmer::from_u64(0, (sshash.k() + 1) as u8);
        sshash.k2u_fw(&km);
    }

    // Test case when all k-mers are in the skew index.
    #[test]
    fn tiny_skew_index() {
        // tiny.seg (k=7)
        // 3	CACACACCAC
        // 2	CCTCAATACG

        let w = 3;
        let skew_param = 0; //
        let unitigs = load_unitigs(TINY_CF_PREFIX);
        let bh = WyHashState::default();
        // All k-mers in skew index
        let skew_index =
            SSHash::from_unitig_set(unitigs.clone(), w, skew_param, bh.clone()).unwrap();

        assert_eq!(skew_index.n_kmers_in_skew_index(), skew_index.n_kmers());
        // No k-mers in skew index
        let no_skew_index = SSHash::from_unitig_set_no_skew_index(unitigs, w, bh).unwrap();

        assert!(skew_index.num_bits() > no_skew_index.num_bits());

        let km1 = CanonicalKmer::from("CACACAC");
        let km2 = CanonicalKmer::from("ACACCAC");
        let km3 = CanonicalKmer::from("CCTCAAT");
        let km4 = CanonicalKmer::from("CCTCAAT");

        let kms = vec![km1, km2, km3, km4];

        for km in kms {
            let lhs = skew_index.k2u(&km).unwrap();
            let rhs = no_skew_index.k2u(&km).unwrap();
            assert_eq!(lhs, rhs);

            let mut rc = km.clone();
            rc.swap();

            let lhs = skew_index.k2u(&rc).unwrap();
            let rhs = no_skew_index.k2u(&rc).unwrap();
            assert_eq!(lhs, rhs);
        }

        let km = CanonicalKmer::from("AAAAAAA");
        let pos = skew_index.k2u(&km);
        assert_eq!(pos, None);
    }
}
