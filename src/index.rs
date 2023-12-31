use kmers::naive_impl::CanonicalKmer;
use serde::{Deserialize, Serialize};

use crate::{
    kphf::{K2UPos, K2U},
    refseq::{RefSeqCollection, RefSeqIter, RefSeqSlice},
    Orientation,
};

pub mod caching;
pub mod defaults;
pub mod dense_unitig_table;
pub mod piscem_index;
pub mod validate;

// mod sparse_unitig_table;
// mod refseq_iter;

// rexports
pub use dense_unitig_table::{DenseUnitigTable, PiscemUnitigTable};
pub use validate::Validate;
// pub use sparse_unitig_table::SparseUnitigTable;

/// Info for a mapped k-mer locating where and in which orientation
/// a mapped k-mer occurs on a reference
#[derive(Debug, Clone, PartialEq)]
pub struct MappedRefPos {
    pub ref_id: usize,
    pub pos: usize,
    pub o: Orientation,
}

impl MappedRefPos {
    pub fn new(ref_id: usize, pos: usize, o: Orientation) -> Self {
        Self { ref_id, pos, o }
    }

    /// New a MappedRefPos mapping in the forward orientation
    pub fn new_fw(ref_id: usize, pos: usize) -> Self {
        Self::new(ref_id, pos, Orientation::Forward)
    }

    /// New a MappedRefPos mapping in the backwars (reverse complement) direction
    pub fn new_rc(ref_id: usize, pos: usize) -> Self {
        Self::new(ref_id, pos, Orientation::Backward)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModIndex<H, T> {
    pub(super) base: BaseIndex,
    u2pos: T,
    k2u: H,
    refs: RefSeqCollection, // NOTE refseq collection can contain ref lengths but not sequences
}

impl<H, T> ModIndex<H, T> {
    pub fn index_type(&self) -> ModIndexType {
        self.base.index_type.clone()
    }
}

impl<H, T> ModIndex<H, T>
where
    H: K2U,
{
    pub fn k(&self) -> usize {
        self.k2u.k()
    }

    pub fn get_ref_names(&self) -> Vec<String> {
        log::warn!("FIX ME");
        Vec::new()
    }

    pub fn as_refseqs(&self) -> &RefSeqCollection {
        &self.refs
    }

    pub fn from_parts(base: BaseIndex, k2u: H, u2pos: T, refs: RefSeqCollection) -> Self {
        Self {
            base,
            u2pos,
            k2u,
            refs,
        }
    }

    pub fn as_u2pos(&self) -> &T {
        &self.u2pos
    }

    pub fn as_k2u(&self) -> &H {
        &self.k2u
    }

    // TODO merge this code block with impl below
    pub fn n_kmers(&self) -> usize {
        self.k2u.n_kmers()
    }

    pub fn n_unitigs(&self) -> usize {
        self.k2u.n_unitigs()
    }

    pub fn unitig_len(&self, uid: usize) -> usize {
        self.k2u.unitig_len(uid)
    }

    pub fn sum_unitigs_len(&self) -> usize {
        self.k2u.sum_unitigs_len()
    }
    pub fn has_refseq(&self) -> bool {
        self.refs.has_seq()
    }

    pub fn n_refs(&self) -> usize {
        self.refs.n_refs()
    }

    pub fn iter_refs(&self) -> RefSeqIter {
        self.refs.iter_refs()
    }
}

#[derive(Clone, PartialEq, Debug)]
pub struct ProjectedHits<OccsT> {
    pub kmer: CanonicalKmer,        // queried kmer
    pub k2upos: K2UPos,             // mapped offset orientation, unitig id, unitig len
    pub encoded_unitig_occs: OccsT, // encoded unitig pos, orientation
}

pub trait GetRefPos {
    type Hits<'a>
    where
        Self: 'a;

    fn get_ref_pos<'a>(&'a self, km: &CanonicalKmer) -> Option<Self::Hits<'a>>;
    fn get_ref_pos_eager(&self, km: &CanonicalKmer) -> Option<Vec<MappedRefPos>> {
        let hits = self.get_ref_pos(km)?;
        Some(self.project_hits(hits))
    }

    fn project_hits<'a>(&'a self, hits: Self::Hits<'a>) -> Vec<MappedRefPos>;
}

impl<H, T> GetRefPos for ModIndex<H, T>
where
    H: K2U,
    T: U2Pos,
{
    type Hits<'a> = ProjectedHits<T::EncodedOccs<'a>>
    where
        Self: 'a;

    fn get_ref_pos<'a>(&'a self, km: &CanonicalKmer) -> Option<Self::Hits<'a>> {
        if km.len() != self.k() {
            panic!(
                "Got query k-mer size k={}, expected k={}.",
                km.len(),
                self.k()
            );
        }

        let k2upos = self.k2u.k2u(km)?;
        let encoded_unitig_occs = self.u2pos.encoded_unitig_occs(k2upos.unitig_id);
        Some(ProjectedHits {
            kmer: km.clone(),
            k2upos,
            encoded_unitig_occs,
        })
    }

    fn project_hits<'a>(&'a self, hits: Self::Hits<'a>) -> Vec<MappedRefPos> {
        let u_occs = self
            .u2pos
            .decode_unitig_occs(&self.k2u, hits.encoded_unitig_occs);
        project_onto_u_occs(self.k(), hits.k2upos, &u_occs)
    }
}

/// project result of K2U query (offeset of k-mer on unitig) onto a list of unitig occurrences
/// retrieving a k-mer's mapped reference positions
pub fn project_onto_u_occs(k: usize, k2upos: K2UPos, u_occs: &[UnitigOcc]) -> Vec<MappedRefPos> {
    u_occs
        .iter()
        .map(|occ| project_onto_u_occ(k, &k2upos, occ))
        .collect()
}

/// project result of K2U query (offeset of k-mer on unitig) onto one unitig occurrence
/// retrieving a k-mer's mapped reference position
#[inline]
pub fn project_onto_u_occ(k: usize, k2upos: &K2UPos, occ: &UnitigOcc) -> MappedRefPos {
    let ref_id = occ.ref_id;
    // let k = self.k();
    let contig_pos = occ.pos;

    let pos = match occ.o {
        Orientation::Forward => k2upos.pos + contig_pos,
        Orientation::Backward => contig_pos + (k2upos.unitig_len - k2upos.pos) - k,
    };

    let o = match k2upos.o {
        kmers::naive_impl::MatchType::IdentityMatch => Orientation::Forward,
        kmers::naive_impl::MatchType::TwinMatch => Orientation::Backward,
        _ => unreachable!(),
    };

    let o = match occ.o {
        Orientation::Forward => o,
        Orientation::Backward => o.reverse(),
    };

    MappedRefPos { ref_id, pos, o }
}

// ----------------------------------------------------------------------------
// BaseIndex

#[derive(Clone, Debug, Serialize, Deserialize, Eq, PartialEq)]
pub struct BaseIndex {
    // basic provenance information
    mazu_version: String,
    index_type: ModIndexType,
    // pub reference_gfa: Vec<String>, //TODO figure out what to do about reference / reduced GFA and provenance info

    // pub sampling_type: PufferfishType,
    // pub k: usize,
    // pub num_kmers: usize,
    // pub num_contigs: usize,
    // pub seq_len: usize,
    // have_ref_seq: bool,
    // pub have_edge_vec: bool,
    // Traversal
    // pub last_seq_pos: usize,
    metadata: Option<IndexMetadata>,
}

impl BaseIndex {
    pub fn new() -> Self {
        Self {
            mazu_version: crate::get_mazu_version(),
            index_type: ModIndexType::unknown(),
            metadata: None,
        }
    }

    pub fn set_index_type(mut self, index_type: ModIndexType) -> Self {
        self.index_type = index_type;
        self
    }

    pub fn set_metadata(mut self, metadata: IndexMetadata) -> Self {
        self.metadata = Some(metadata);
        self
    }
}

impl Default for BaseIndex {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Eq, PartialEq)]
pub struct IndexMetadata {
    pub have_edge_vec: bool,
    pub seq_hash: String,
    pub name_hash: String,
    pub seq_hash_512: String,
    pub name_hash_512: String,
    pub decoy_seq_hash: String,
    pub decoy_name_hash: String,
    pub num_decoys: usize,
    pub first_decoy_index: usize,
    pub keep_duplicates: bool,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub enum ModIndexType {
    PF1Dense,  // Dense index deserialized from pufferfish (C++)
    PF1Sparse, // Sparse index deserialized from pufferfish (C++)
    PufferfishDense,
    // PufferfishSparse,
    Piscem,
    // SparseRS(PosSamplingConfig),
    // SparseContigDensePos(ContigSamplingStrategy),
    // SparseContigSparsePos(ContigSamplingStrategy),
    // // SparseContigSparsePos(ContigSamplingStrategy, PosSamplingConfig),
    // SparseDenseV2(ContigSamplingStrategy),
    // SparseSparseV2(ContigSamplingStrategy),
    Custom(String), // Custom index type with description
}

impl ModIndexType {
    fn unknown() -> Self {
        Self::Custom("Custom, unspecified, modular index type".to_string())
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnitigOcc {
    pub ref_id: usize,
    pub pos: usize,
    pub o: Orientation,
}

impl UnitigOcc {
    #[inline]
    pub fn is_forward(&self) -> bool {
        matches!(self.o, Orientation::Forward)
    }

    /// Encode [UnitigOcc] into a [u64]word (as in pufferfish),
    ///
    /// pufferfish encoding:
    // | 1 bit orientation | 31 bit position | 32 bit ref ID |
    #[inline]
    pub fn encode_pf1(&self) -> u64 {
        debug_assert!(self.ref_id < u32::MAX as usize);
        debug_assert!(self.pos < 0x80000000);
        let mut word = self.pos as u64;
        if self.is_forward() {
            word |= 0x80000000
        };

        word <<= 32;
        word |= self.ref_id as u64;
        word
    }

    /// Decode [UnitigOcc] from a [u64] word (as in pufferfish),
    #[inline]
    pub fn decode_pf1(word: u64) -> Self {
        let ref_id = (word & 0xFFFFFFFF) as usize;
        let pos_word = word >> 32;
        let pos = (pos_word & 0x7FFFFFFF) as usize;
        let o = if (pos_word & 0x80000000) > 0 {
            Orientation::Forward
        } else {
            Orientation::Backward
        };

        Self { ref_id, pos, o }
    }
}

pub trait U2Pos {
    // trait bounds for Clone and Debug to make GAT make sense,
    // otherwise U2Pos becomes clumsy to use, requiring bounds for
    // Clone and Debug on the GAT
    type EncodedOccs<'a>: Clone + std::fmt::Debug
    where
        Self: 'a;
    fn encoded_unitig_occs(&self, ui: usize) -> Self::EncodedOccs<'_>;
    fn decode_unitig_occs<T>(&self, meta: &T, occs: Self::EncodedOccs<'_>) -> Vec<UnitigOcc>;
    fn n_unitig_occs(&self, occs: Self::EncodedOccs<'_>) -> usize;

    fn empty_encoded_occs(&self) -> Self::EncodedOccs<'_>;
}

impl<H: K2U, T> ModIndex<H, T> {
    pub fn iter_unitigs_on_ref<'a, 'b>(
        &'a self,
        refseq: &RefSeqSlice<'b>,
    ) -> RefSeqContigIterator<'a, 'b, H> {
        let end_pos = refseq.len() - self.k() + 1;

        RefSeqContigIterator {
            pos: 0,
            refseq: refseq.clone(),
            end_pos,
            k2u: &self.k2u,
        }
    }
}

#[derive(Debug)]
pub struct RefSeqContigIterator<'a, 'b, K> {
    // start position of current contig
    pos: usize,
    refseq: RefSeqSlice<'b>,
    end_pos: usize,
    k2u: &'a K,
}

#[derive(Debug)]
pub struct RefSeqUnitigOcc {
    pub unitig_id: usize,
    pub unitig_len: usize,
    pub pos: usize,
    pub o: Orientation,
}

impl<'a, 'b, K: K2U> Iterator for RefSeqContigIterator<'a, 'b, K> {
    type Item = RefSeqUnitigOcc;
    fn next(&mut self) -> Option<Self::Item> {
        if self.pos < self.end_pos {
            let km = self.refseq.get_kmer(self.k2u.k(), self.pos);
            let km = CanonicalKmer::from(km);
            let hit = self.k2u.k2u(&km).unwrap();
            let o = match hit.o {
                kmers::naive_impl::MatchType::IdentityMatch => Orientation::Forward,
                kmers::naive_impl::MatchType::TwinMatch => Orientation::Backward,
                _ => unreachable!(),
            };

            let occ = RefSeqUnitigOcc {
                // ref_id: self.refseq.ref_id(),
                unitig_id: hit.unitig_id,
                unitig_len: hit.unitig_len,
                pos: self.pos,
                o,
            };

            self.pos += occ.unitig_len - self.k2u.k() + 1;

            Some(occ)
        } else {
            None
        }
    }
}
