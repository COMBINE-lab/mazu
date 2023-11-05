use kmers::naive_impl::{
    seq_vector::{SeqVector, SeqVectorSlice},
    Kmer,
};
use serde::{Deserialize, Serialize};
use simple_sds::{
    bit_vector::BitVector,
    ops::{BitVec, Rank, Select, PredSucc},
    sparse_vector::{SparseBuilder, SparseVector},
    raw_vector::{AccessRaw, RawVector},
};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use crate::{
    cuttlefish::{CfFiles, CfInfo},
    elias_fano::EFVector,
    Result,
};


/// Compact encoding for a collection of unitigs 
/// but without the unitig sequence itself.
///
/// This contains the unitig lengths and (indexed) offsets as well 
/// as the value of k for which the unitigs were built.
///
/// This is essentially a `UnitigSet`, but without the sequence 
/// itself.
#[derive(Debug)]
pub struct UnitigSetInfo {
    pub(crate) k : usize,
    pub(crate) accum_lens : SparseVector,
}

impl UnitigSetInfo {
    /// Create a `UnitigSetInfo` struct from `cuttlefish` output
    ///
    /// Consumes the files in "reduced GFA format" from `cuttlefish`.
    /// Returns a UnitigSet, and a HashMap mapping IDs given by `cuttlefish`
    /// (that are not sequential) to sequention IDs [0..N] for N unitig sequences.
    pub fn from_cf_reduced_gfa(cf_files: &CfFiles) -> Result<(Self, HashMap<usize, usize>)> {
        let cf_info = CfInfo::from_path(&cf_files.json)?;
        let len = cf_info.total_len();

        // 1) open cf segs and insert the segments
        let mut uid_to_idx = HashMap::new();
        // the +1 below is because to ease working with the 
        // resulting bitvector, we will place a 1 at the end 
        // of the last unitig as well.
        let mut accum_lens = SparseBuilder::new(
            len + 1, cf_info.n_unitigs())
            .expect("construct accum_lens builder");

        let f = File::open(&cf_files.segs).unwrap();
        let f = BufReader::new(f);

        let mut prefix_sum = 0;

        for (i, line) in f.lines().enumerate() {
            let line = line?;
            let (id, seq) = line.split_once('\t').expect("Cannot split .cf_seqs line");
            let id: usize = id.parse().expect("Failed to parse unitig ID");
            let ulen = seq.len();

            uid_to_idx.insert(id, i);
            prefix_sum += ulen;
            accum_lens.set(prefix_sum - 1);
        }

        let accum_lens = SparseVector::try_from(accum_lens)
            .expect("building SparseVector accum_lens from builder.");

        Ok((
            Self {
                k: cf_info.k(),
                accum_lens,
            },
            uid_to_idx,
        ))
    }

    /// Map the given _global_ position on concatenated unitig sequences to
    /// the corresponding unitig ID. This leaky abstraction is useful for
    /// k-mer indexing. (See [crate::kphf::SSHash] and [crate::kphf::PFHash])
   #[inline]
    pub fn pos_to_id(&self, pos: usize) -> usize {
        self.accum_lens.rank(pos)
    }

    /// Number of unitigs encoded in this set
    pub fn n_unitigs(&self) -> usize {
        self.accum_lens.count_ones()
    }

    /// Return length of unitig `i`
    pub fn unitig_len(&self, i: usize) -> usize {
        if i > 0 {
            (self.accum_lens.select(i).unwrap() - self.accum_lens.select(i-1).unwrap()) as usize
        } else {
            self.accum_lens.select(0).unwrap() + 1
        }
    }

    /// Get where sequence unitig `i` starts on global cancatenated [SeqVector]
    pub fn unitig_start_pos(&self, i: usize) -> usize {
        if i > 0 {
            self.accum_lens.select(i - 1).unwrap() + 1 as usize
        } else {
            0_usize
        }
    }

    /// Get where sequence unitig `i` ends on global cancatenated [SeqVector]
    /// This function returns 1 index beyond the last position of the unitig 
    /// (that is, we employ half-open range indexing so that unitig i spans
    /// [`unitig_start_pos(i)`, `unitig_end_pos(i)`) )
    pub fn unitig_end_pos(&self, i: usize) -> usize {
        self.accum_lens.select(i).unwrap() + 1 as usize
    }

    /// Return the total length of the [UnitigSetInfo] ()
    pub fn total_len(&self) -> usize {
        self.accum_lens.select(self.n_unitigs() - 1).unwrap() + 1 as usize
    }

    /// Return number of k-mers encoded in [UnitigSetInfo]
    pub fn n_kmers(&self) -> usize {
        self.total_len() - (self.k * self.n_unitigs()) + self.n_unitigs()
    }

    pub fn k(&self) -> usize {
        self.k
    }
}


/// Compact encoding for a collection of unitigs
///
/// [UnitigSet] concatenates unitig sequences into a 2-bit encoded
/// [sequence vector](SeqVector). Unitig-Unitig boundaries are
/// demarcated by a rank/select supported bitvector.
// NOTE: update Self::num_bits if struct changes
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct UnitigSet {
    pub(crate) k: usize,
    pub(crate) useq: SeqVector,
    pub(crate) accum_lens: EFVector,
    pub(crate) bv: BitVector,
}

impl UnitigSet {
    fn accum_len_bits(&self) -> usize {
        // self.accum_lens.len() * std::mem::size_of::<usize>() * 8 // use len Vec::len() instead of capacity
        self.accum_lens.num_bits()
    }
    pub fn num_bits(&self) -> usize {
        std::mem::size_of::<usize>() * 8
            + self.useq.num_bits()
            + self.accum_len_bits()
            + self.bv.num_bits()
    }

    pub fn print_stats(&self) {
        log::info!(
            "{:.5} bits per kmer",
            self.num_bits() as f64 / self.n_kmers() as f64
        );
        log::info!(
            "seq: {:.5} bits per nuc",
            self.useq.num_bits() as f64 / self.total_len() as f64
        );
        log::info!(
            "bv: {:.5} bits per nuc",
            self.bv.num_bits() as f64 / self.total_len() as f64
        );
        log::info!(
            "accum_lens: {:.5} bits per km",
            self.accum_len_bits() as f64 / self.n_kmers() as f64
        );
        log::info!(
            "            {:.5} bits per nuc",
            self.accum_len_bits() as f64 / self.total_len() as f64
        );
    }

    /// Create a `UnitigSet` from a list of strings
    pub fn from_seqs(seqs: &[String], k: usize) -> crate::err::Result<Self> {
        let len = seqs.iter().map(|s| s.len()).sum();

        let mut prefix_sum = 0;
        let mut useq = SeqVector::with_capacity(len);
        let mut accum_lens = Vec::with_capacity(seqs.len() + 1);

        for seq in seqs {
            let ulen = seq.len();

            useq.push_chars(seq.as_bytes());
            accum_lens.push(prefix_sum);
            prefix_sum += ulen;
        }
        accum_lens.push(prefix_sum);

        let mut bv = RawVector::with_len(prefix_sum, false);
        for l in &accum_lens[1..] {
            bv.set_bit(l - 1, true);
        }

        let mut bv = BitVector::from(bv);
        bv.enable_rank();

        let accum_lens = EFVector::from_usize_slice(&accum_lens)?;

        Ok(Self {
            k,
            useq,
            accum_lens,
            bv,
        })
    }

    /// Create a `UnitigSet` from `cuttlefish` prefix
    pub fn from_cf_prefix<P: AsRef<Path>>(prefix: P) -> Result<(Self, HashMap<usize, usize>)> {
        let cf_files = CfFiles::new(prefix);
        Self::from_cf_reduced_gfa(&cf_files)
    }

    /// Create a `UnitigSet` from `cuttlefish` output
    ///
    /// Consumes the files in "reduced GFA format" from `cuttlefish`.
    /// Returns a UnitigSet, and a HashMap mapping IDs given by `cuttlefish`
    /// (that are not sequential) to sequention IDs [0..N] for N unitig sequences.
    pub fn from_cf_reduced_gfa(cf_files: &CfFiles) -> Result<(Self, HashMap<usize, usize>)> {
        let cf_info = CfInfo::from_path(&cf_files.json)?;

        let len = cf_info.total_len();
        // 1) open cf segs and insert the segments
        let mut useq = SeqVector::with_capacity(len);
        let mut uid_to_idx = HashMap::new();
        let mut accum_lens = Vec::with_capacity(cf_info.n_unitigs() + 1);

        let f = File::open(&cf_files.segs).unwrap();
        let f = BufReader::new(f);

        let mut prefix_sum = 0;

        for (i, line) in f.lines().enumerate() {
            let line = line?;
            let (id, seq) = line.split_once('\t').expect("Cannot split .cf_seqs line");
            let id: usize = id.parse().expect("Failed to parse unitig ID");
            let ulen = seq.len();

            useq.push_chars(seq.as_bytes());
            uid_to_idx.insert(id, i);
            accum_lens.push(prefix_sum);
            prefix_sum += ulen;
        }
        accum_lens.push(prefix_sum);

        let mut bv = RawVector::with_len(prefix_sum, false);
        for l in &accum_lens[1..] {
            bv.set_bit(l - 1, true);
        }

        let mut bv = BitVector::from(bv);
        bv.enable_rank();

        let accum_lens = EFVector::from_usize_slice(&accum_lens)?;

        Ok((
            Self {
                k: cf_info.k(),
                useq,
                accum_lens,
                bv,
            },
            uid_to_idx,
        ))
    }

    /// Return the encoded number of unitigs
    pub fn n_unitigs(&self) -> usize {
        self.accum_lens.len() - 1
    }

    /// Returns a iterator over canonical k-mer iterators for each unitig sequence
    pub fn chunked_unitig_canonical_kmers(&self) -> ChunkedUnitigSetCanonicalKmers {
        ChunkedUnitigSetCanonicalKmers { unitigs: self }
    }

    /// Return length of unitig `i`
    pub fn unitig_len(&self, i: usize) -> usize {
        (self.accum_lens.get(i + 1) - self.accum_lens.get(i)) as usize
    }

    /// Map the given _global_ position on concatenated unitig sequences to
    /// the corresponding unitig ID. This leaky abstraction is useful for
    /// k-mer indexing. (See [crate::kphf::SSHash] and [crate::kphf::PFHash])
    pub fn pos_to_id(&self, pos: usize) -> usize {
        self.bv.rank(pos)
    }

    /// Get [SeqVectorSlice] corresponding to unitig `i`.
    pub fn unitig_seq(&self, i: usize) -> SeqVectorSlice {
        let start = self.accum_lens.get(i) as usize;
        let end = self.accum_lens.get(i + 1) as usize;
        self.useq.slice(start, end)
    }

    /// Get where sequence unitig `i` starts on global cancatenated [SeqVector]
    pub fn unitig_start_pos(&self, i: usize) -> usize {
        self.accum_lens.get(i) as usize
    }

    /// Get where sequence unitig `i` ends on global cancatenated [SeqVector]
    pub fn unitig_end_pos(&self, i: usize) -> usize {
        self.accum_lens.get(i + 1) as usize
    }

    /// Return the total length of the [UnitigSet] ()
    pub fn total_len(&self) -> usize {
        self.accum_lens.get(self.n_unitigs()) as usize
    }

    /// Return number of k-mers encoded in [UnitigSet]
    pub fn n_kmers(&self) -> usize {
        self.total_len() - (self.k * self.n_unitigs()) + self.n_unitigs()
    }

    /// Get k-mer from given position on _concatenated_ unitig sequences.
    ///
    /// This leaky abstraction is useful for k-mer indexing.
    /// (See [crate::kphf::SSHash] and [crate::kphf::PFHash])
    pub fn get_kmer_from_useq_pos(&self, pos: usize) -> Kmer {
        let km = self.get_kmer_u64_from_useq_pos(pos);
        Kmer::from_u64(km, self.k as u8)
    }

    /// Get k-mer (as `u64` word) from given position on _concatenated_ unitig sequences.
    pub fn get_kmer_u64_from_useq_pos(&self, pos: usize) -> u64 {
        // May return an invalid kmer on unitig-unitig boundary on useq
        self.useq.get_kmer_u64(pos, self.k)
    }

    /// Checks if a position on _concatenated_ unitig sequence is a position of a valid k-mer.
    ///
    /// I.e., returns false if a k-mer at given position lies on unitig-unitig boundary on
    /// concatenated unitig-sequence.
    pub fn is_valid_useq_pos(&self, pos: usize) -> bool {
        let last_kmer_pos = self.total_len() - self.k();
        if pos > last_kmer_pos {
            false
        } else {
            // check if kmer at pos crosses unitig-unitig boundary
            // the *ends* of unitigs are marked with 1. so the no valid kmer has k-1 prefix with a set bit in bv
            let word = unsafe { self.bv.as_ref().int(pos, self.k() - 1) };
            word == 0
        }
    }

    pub fn k(&self) -> usize {
        self.k
    }
}

// Iterators for MPHF construction
// Chunked iterators iterate over unitig iterators generating canonicalized k-mer seqs as u64 words

/// Iterator over [kmers::naive_impl::CanonicalKmer]s of a [SeqVectorSlice]
#[derive(Clone, Debug)]
pub struct SeqVecCanonicalKmerIterU64<'a> {
    sv: SeqVectorSlice<'a>,
    k: usize,
    curr: usize,
}

/// Iterator over [kmers::naive_impl::CanonicalKmer]s of a [SeqVectorSlice]
#[derive(Clone, Debug)]
pub struct SeqVecSliceCanonicalKmers<'a> {
    sv: SeqVectorSlice<'a>,
    k: usize,
}

/// An `IntoIterator` for the chunked iterator that iterates over canonical k-mer iterators
/// over each unitig sequence
#[derive(Debug, Clone)]
pub struct ChunkedUnitigSetCanonicalKmers<'a> {
    unitigs: &'a UnitigSet,
}

/// A chunked iterator that iterates over canonical k-mer iterators over each unitig sequence
#[derive(Debug, Clone)]
pub struct ChunkedUnitigIterator<'a> {
    unitigs: &'a UnitigSet,
    curr: usize,
}

impl<'a> IntoIterator for SeqVecSliceCanonicalKmers<'a> {
    type Item = u64;
    type IntoIter = SeqVecCanonicalKmerIterU64<'a>;

    fn into_iter(self) -> Self::IntoIter {
        Self::IntoIter {
            sv: self.sv,
            k: self.k,
            curr: 0,
        }
    }
}

impl<'a> Iterator for SeqVecCanonicalKmerIterU64<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr < (self.sv.len() - self.k + 1) {
            let km = self.sv.get_kmer(self.curr, self.k);
            let km = km.to_canonical();
            let res = Some(km.into_u64());
            self.curr += 1;
            res
        } else {
            None
        }
    }
}

impl ExactSizeIterator for SeqVecCanonicalKmerIterU64<'_> {
    fn len(&self) -> usize {
        self.sv.len() - self.k + 1 - self.curr
    }
}

impl<'a> IntoIterator for &ChunkedUnitigSetCanonicalKmers<'a> {
    type Item = SeqVecSliceCanonicalKmers<'a>;
    type IntoIter = ChunkedUnitigIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        Self::IntoIter {
            unitigs: self.unitigs,
            curr: 0,
        }
    }
}

impl<'a> Iterator for ChunkedUnitigIterator<'a> {
    type Item = SeqVecSliceCanonicalKmers<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.curr < self.unitigs.n_unitigs() {
            let slice = self.unitigs.unitig_seq(self.curr);
            let iter = Self::Item {
                sv: slice,
                k: self.unitigs.k(),
            };
            let res = Some(iter);
            self.curr += 1;
            res
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {
    use crate::test_utils::*;
    use crate::unitig_set::UnitigSetInfo;
    use crate::unitig_set::CfFiles;

    #[test]
    fn unitigs() {
        // Unitig sequences in tiny.seq
        // 3	CACACACCAC
        // 2	CCTCAATACG

        let unitigs = load_unitigs(TINY_CF_PREFIX);

        assert_eq!(unitigs.k(), 7);

        let u0 = unitigs.unitig_seq(0);
        assert_eq!("CACACACCAC", u0.to_string());

        let u1 = unitigs.unitig_seq(1);
        assert_eq!("CCTCAATACG", u1.to_string());

        assert_eq!(unitigs.unitig_len(0), 10);
        assert_eq!(unitigs.unitig_len(1), 10);

        for i in 0..10 {
            assert_eq!(unitigs.pos_to_id(i), 0);
        }

        for i in 10..20 {
            assert_eq!(unitigs.pos_to_id(i), 1);
        }

        assert_eq!(unitigs.total_len(), 20);
    }


    use simple_sds::ops::{BitVec, Rank, Select, SelectZero, PredSucc};

    #[test]
    fn unitig_info() {
        let cf_files = CfFiles::new(TINY_CF_PREFIX);
        if let Ok((unitig_set_info, _cfid_2_uid)) = UnitigSetInfo::from_cf_reduced_gfa(&cf_files) {
            assert_eq!(unitig_set_info.k(), 7);

            assert_eq!(unitig_set_info.unitig_len(0), 10);
            assert_eq!(unitig_set_info.unitig_len(1), 10);

            assert_eq!(unitig_set_info.unitig_start_pos(0), 0);
            assert_eq!(unitig_set_info.unitig_start_pos(1), 10);

            assert_eq!(unitig_set_info.unitig_end_pos(0), 10);
            assert_eq!(unitig_set_info.unitig_end_pos(1), 20);

            assert_eq!(unitig_set_info.total_len(), 20);

            for i in 0..10 {
                assert_eq!(unitig_set_info.pos_to_id(i), 0);
            }

            for i in 10..20 {
                assert_eq!(unitig_set_info.pos_to_id(i), 1);
            }

        } else {
            panic!("couldn't construct the UnitigSetInfo on the tiny example data.");
        }
    }


    #[test]
    fn unitig_parsers_agree() {
        let unitigs = load_unitigs(YEAST_CF_PREFIX);
        let cf_files = CfFiles::new(YEAST_CF_PREFIX);
        if let Ok((unitig_set_info, _cfid_2_uid)) = UnitigSetInfo::from_cf_reduced_gfa(&cf_files) {
            assert_eq!(unitig_set_info.k(), unitigs.k());

            let n_kmers = unitigs.n_kmers();
            assert_eq!(unitig_set_info.n_kmers(), n_kmers);
 
            let n_unitigs = unitigs.n_unitigs();
            assert_eq!(unitig_set_info.n_unitigs(), n_unitigs);

            for i in 0..n_unitigs {
                assert_eq!(unitig_set_info.unitig_len(i), unitigs.unitig_len(i));
                assert_eq!(unitig_set_info.unitig_start_pos(i), unitigs.unitig_start_pos(i));
                assert_eq!(unitig_set_info.unitig_end_pos(i), unitigs.unitig_end_pos(i),
                    "testing the end position of unitig {}", i);
            }

            let tlen = unitigs.total_len();
            assert_eq!(unitig_set_info.total_len(), tlen);

            for i in 0..tlen {
                assert_eq!(unitig_set_info.pos_to_id(i), unitigs.pos_to_id(i));
            }
        } else {
            panic!("couldn't construct the UnitigSetInfo on the tiny example data.");
        }
    }
}
