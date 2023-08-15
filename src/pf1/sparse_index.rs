use kmers::naive_impl::seq_vector::SeqVector;

use simple_sds::{int_vector::IntVector, ops::BitVec};

use super::{boophf::BooPHF, cpp::*, PF1Info};
use crate::{
    elias_fano::EFVector,
    index::{BaseIndex, DenseUnitigTable, ModIndex, PufferfishType},
    kphf::SampledPFHash,
    pf1::PF1FilePaths,
    refseq::RefSeqCollection,
    unitig_set::UnitigSet,
    Error, Result,
};
use log::debug;
use simple_sds::{
    bit_vector::BitVector,
    ops::{Rank, Select},
};
use std::path::Path;

#[allow(non_camel_case_types)]
type pf1_hash_t = SampledPFHash<BooPHF<u64>>;
pub type SparseIndex = ModIndex<pf1_hash_t, DenseUnitigTable>;

impl AsRef<UnitigSet> for SparseIndex {
    fn as_ref(&self) -> &UnitigSet {
        self.as_k2u().as_ref()
    }
}

impl DeserializeFromCpp for SparseIndex {
    fn deserialize_from_cpp<P: AsRef<Path>>(dir: P) -> Result<Self> {
        debug!("Loading base index");
        if let Ok(base) = BaseIndex::deserialize_from_cpp(&dir) {
            let files = PF1FilePaths::new(&dir);

            debug!("Loading Kmer PHF");
            debug!("* Loading seq");
            let useq = SeqVector::from_compact_serialized(files.seq)?;

            debug!("* Loading bv");
            let bv = {
                let mut bv = BitVector::from_compact_serialized(files.rank)?;
                bv.enable_rank();
                bv.enable_select();
                bv
            };
            debug!("* Constructing unitig set");

            log::warn!("* Computing prefix sums for accumulated lenghts, need to deprecate this");

            let n_unitigs = bv.count_ones();
            let mut accum_lens = Vec::with_capacity(n_unitigs);
            let accum = 0;
            accum_lens.push(accum);
            for ui in 1..(n_unitigs + 1) {
                let start_pos = bv.select(ui - 1).unwrap() + 1;
                accum_lens.push(start_pos);
            }

            let accum_lens = EFVector::from_usize_slice(&accum_lens)?;

            let unitigs = UnitigSet {
                k: base.k,
                useq,
                accum_lens,
                bv,
            };

            debug!("* Loading MPHF and other kmer hash data structures");
            let mphf = BooPHF::<u64>::deserialize_from_cpp(files.mphf)?;
            let sampled_pos = IntVector::from_compact_serialized(files.sample_pos)?;
            let canonical_vec = BitVector::from_compact_serialized(files.canonical)?;
            let direction_vec = BitVector::from_compact_serialized(files.direction)?;
            let ext_sizes = IntVector::from_compact_serialized(files.extension_lengths)?;
            let ext_bases = IntVector::from_compact_serialized(files.extension_bases)?;
            let sampled_vec = {
                let mut sv = BitVector::from_compact_serialized(files.presence)?;
                sv.enable_rank();
                sv.enable_select();
                sv
            };

            let info = PF1Info::load(files.info_json);
            let sample_size = info.sample_size.unwrap();
            let extension_size = info.extension_size.unwrap();

            let k2u = SampledPFHash {
                unitigs,
                mphf,
                sampled_vec,
                canonical_vec,
                direction_vec,
                ext_sizes,
                ext_bases,
                sampled_pos,
                sample_size,
                extension_size,
            };

            debug!("Loading contig table");
            let ctg_table = DenseUnitigTable::deserialize_from_cpp(&dir)?;

            let refs = RefSeqCollection::deserialize_from_cpp(&dir)?;
            let u2pos = ctg_table;
            Ok(Self::from_parts(base, k2u, u2pos, refs))
        } else {
            Err(Error::IndexLoad)
        }
    }
}

impl SparseIndex {
    pub fn to_pf1_info(&self) -> PF1Info {
        let base = &self.base;
        PF1Info {
            sampling_type: PufferfishType::SparseRS,
            index_version: base.index_version,
            reference_gfa: base.reference_gfa.clone(),
            kmer_size: base.k,
            num_kmers: self.n_kmers(),
            num_contigs: self.n_unitigs(),
            seq_len: self.sum_unitigs_len(),
            have_ref_seq: self.has_refseq(),
            have_edge_vec: base.have_edge_vec,
            seq_hash: base.seq_hash.clone(),
            name_hash: base.name_hash.clone(),
            seq_hash_512: base.seq_hash_512.clone(),
            name_hash_512: base.name_hash_512.clone(),
            decoy_seq_hash: base.decoy_seq_hash.clone(),
            decoy_name_hash: base.decoy_name_hash.clone(),
            num_decoys: base.num_decoys,
            sample_size: Some(self.as_k2u().sample_size), // Putting into enum would break c++ compatiblity)
            extension_size: Some(self.as_k2u().extension_size),
            first_decoy_index: base.first_decoy_index,
            keep_duplicates: base.keep_duplicates,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        kphf::{SSHash, WyHashState, K2U},
        pf1::test_utils::*,
        Validate,
    };

    #[test]
    fn pufferfish_type() {
        let p = to_abs_path(SMALL_TXOME_SPARSE_INDEX);
        let pi = SparseIndex::deserialize_from_cpp(p).unwrap();
        assert_eq!(pi.to_pf1_info().sampling_type, PufferfishType::SparseRS);
        assert_eq!(pi.to_pf1_info().sample_size.unwrap(), 9);
        assert_eq!(pi.to_pf1_info().extension_size.unwrap(), 4);
    }

    #[test]
    fn k2u_tiny() {
        let p = to_abs_path(SMALL_TXOME_SPARSE_INDEX);
        let pi = SparseIndex::deserialize_from_cpp(p).unwrap();
        pi.as_k2u().validate_self()
    }

    #[test]
    fn tiny() {
        let p = to_abs_path(SMALL_TXOME_SPARSE_INDEX);
        let pi = SparseIndex::deserialize_from_cpp(p).unwrap();
        pi.validate_self()
    }
    #[test]
    fn sshash_drop_in() {
        let p = to_abs_path(SMALL_TXOME_SPARSE_INDEX);

        let pi = SparseIndex::deserialize_from_cpp(p).unwrap();
        let unitig_set = pi.as_ref().clone();
        let sshash =
            SSHash::from_unitig_set_no_skew_index(unitig_set, 2, WyHashState::default()).unwrap();

        let pi = ModIndex::from_parts(
            pi.base.clone(),
            sshash,
            pi.as_u2pos().clone(),
            pi.as_refseqs().clone(),
        );
        pi.validate_self()
    }
}
