use kmers::naive_impl::seq_vector::SeqVector;

use simple_sds::{int_vector::IntVector, ops::BitVec};

use log::debug;
use simple_sds::{
    bit_vector::BitVector,
    ops::{Rank, Select, Vector},
};
use std::path::Path;

use super::{boophf::BooPHF, cpp::*, PF1Info};
use crate::{
    elias_fano::EFVector,
    index::{BaseIndex, DenseUnitigTable, ModIndex, PufferfishType},
    kphf::PFHash,
    pf1::PF1FilePaths,
    refseq::RefSeqCollection,
    unitig_set::UnitigSet,
    Error, Result,
};

#[allow(non_camel_case_types)]
type pf1_hash_t = PFHash<BooPHF<u64>>;
pub type DenseIndex = ModIndex<pf1_hash_t, DenseUnitigTable>;

impl AsRef<UnitigSet> for DenseIndex {
    fn as_ref(&self) -> &UnitigSet {
        self.as_k2u().as_ref()
    }
}

impl DenseIndex {
    pub fn to_pf1_info(&self) -> PF1Info {
        let base = &self.base;
        PF1Info {
            sampling_type: PufferfishType::DenseRS,
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
            sample_size: None, // Putting into enum would break c++ compatiblity)
            extension_size: None,
            first_decoy_index: base.first_decoy_index,
            keep_duplicates: base.keep_duplicates,
        }
    }
}

impl DeserializeFromCpp for DenseIndex {
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

            // dbg!(accum);
            // dbg!(n_unitigs);
            let accum_lens = EFVector::from_usize_slice(&accum_lens)?;

            let unitigs = UnitigSet {
                k: base.k,
                useq,
                accum_lens,
                bv,
            };

            debug!("* Loading mphf");
            let mphf = BooPHF::<u64>::deserialize_from_cpp(files.mphf)?;

            debug!("* Loading pos");
            let pos = IntVector::from_compact_serialized(files.pos)?;
            assert_eq!(pos.len(), unitigs.n_kmers());

            let k2u = PFHash::from_parts(unitigs, mphf, pos);

            // debug!("Loading seq");
            // let seq = SeqVector::from_compact_serialized(files.seq)?;
            // assert_eq!(seq.len(), info.seq_len);
            // let kphf = PFHash

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

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        index::GetRefPos,
        kphf::{SSHash, WyHashState},
        pf1::test_utils::*,
        MappedRefPos, Validate,
    };
    use kmers::naive_impl::CanonicalKmer;

    #[test]
    fn pufferfish_type() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(p).unwrap();
        assert_eq!(pi.to_pf1_info().sampling_type, PufferfishType::DenseRS);
    }

    #[test]
    fn tiny_n_unitigs() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(p).unwrap();
        assert_eq!(pi.n_unitigs(), 1);
    }

    #[test]
    fn yeast_n_unitigs() {
        let p = to_abs_path(YEAST_CHR01_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(p).unwrap();
        assert_eq!(pi.n_unitigs(), 577);
    }

    #[test]
    fn tiny_query_does_not_exist() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();

        let kmers = vec!["tat", "ata", "act", "ctg", "cct"];
        let kmers = kmers.iter().map(|s| CanonicalKmer::from(*s));
        // let kmers = kmers.iter().map(|s| Kmer::from(*s));

        for k in kmers {
            assert_eq!(pi.get_ref_pos(&k), None);
        }
    }

    /******************************************************************************/
    // Tiny example: invalid queries panic
    /******************************************************************************/
    #[test]
    #[should_panic]
    fn tiny_kmer_too_big() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(p).unwrap();
        pi.get_ref_pos(&CanonicalKmer::from("aaaaaa"));
    }

    // TODO test panics for k2u queries
    #[test]
    #[should_panic]
    fn tiny_kmer_too_small() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(p).unwrap();
        pi.get_ref_pos(&CanonicalKmer::from("aa"));
    }

    /******************************************************************************/
    // Tiny example: check kmer positions on seq
    /******************************************************************************/
    #[test]
    fn test_tiny_kmer_positions() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(p).unwrap();

        // Indexed string: // AAACCC
        let aaa_pos = vec![MappedRefPos::new_fw(0, 0)];
        let aac_pos = vec![MappedRefPos::new_fw(0, 1)];
        let acc_pos = vec![MappedRefPos::new_fw(0, 2)];
        let ccc_pos = vec![MappedRefPos::new_fw(0, 3)];

        // binary representations;
        let mut aaa = CanonicalKmer::from("aaa");
        let mut aac = CanonicalKmer::from("aac");
        let mut acc = CanonicalKmer::from("acc");
        let mut ccc = CanonicalKmer::from("ccc");

        assert_eq!(pi.get_ref_pos_eager(&aaa), Some(aaa_pos));

        assert_eq!(pi.get_ref_pos_eager(&aac), Some(aac_pos));

        assert_eq!(pi.get_ref_pos_eager(&acc), Some(acc_pos));

        assert_eq!(pi.get_ref_pos_eager(&ccc), Some(ccc_pos));

        let aaa_pos = vec![MappedRefPos::new_rc(0, 0)];
        let aac_pos = vec![MappedRefPos::new_rc(0, 1)];
        let acc_pos = vec![MappedRefPos::new_rc(0, 2)];
        let ccc_pos = vec![MappedRefPos::new_rc(0, 3)];
        aaa.swap();
        aac.swap();
        acc.swap();
        ccc.swap();

        assert_eq!(pi.get_ref_pos_eager(&aaa), Some(aaa_pos));

        assert_eq!(pi.get_ref_pos_eager(&aac), Some(aac_pos));

        assert_eq!(pi.get_ref_pos_eager(&acc), Some(acc_pos));

        assert_eq!(pi.get_ref_pos_eager(&ccc), Some(ccc_pos));
    }

    /******************************************************************************/
    // Tiny example: check kmer positions on seq
    /******************************************************************************/
    #[test]
    fn test_tiny_sshash_k2u() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(p).unwrap();
        let unitig_set = pi.as_ref().clone();
        let sshash =
            SSHash::from_unitig_set_no_skew_index(unitig_set, 2, WyHashState::default()).unwrap();

        let pi = ModIndex::from_parts(
            pi.base.clone(),
            sshash,
            pi.as_u2pos().clone(),
            pi.as_refseqs().clone(),
        );

        // Indexed string: // AAACCC
        let aaa_pos = vec![MappedRefPos::new_fw(0, 0)];
        let aac_pos = vec![MappedRefPos::new_fw(0, 1)];
        let acc_pos = vec![MappedRefPos::new_fw(0, 2)];
        let ccc_pos = vec![MappedRefPos::new_fw(0, 3)];

        // binary representations;
        let mut aaa = CanonicalKmer::from("aaa");
        let mut aac = CanonicalKmer::from("aac");
        let mut acc = CanonicalKmer::from("acc");
        let mut ccc = CanonicalKmer::from("ccc");

        assert_eq!(pi.get_ref_pos_eager(&aaa), Some(aaa_pos));

        assert_eq!(pi.get_ref_pos_eager(&aac), Some(aac_pos));

        assert_eq!(pi.get_ref_pos_eager(&acc), Some(acc_pos));

        assert_eq!(pi.get_ref_pos_eager(&ccc), Some(ccc_pos));

        let aaa_pos = vec![MappedRefPos::new_rc(0, 0)];
        let aac_pos = vec![MappedRefPos::new_rc(0, 1)];
        let acc_pos = vec![MappedRefPos::new_rc(0, 2)];
        let ccc_pos = vec![MappedRefPos::new_rc(0, 3)];
        aaa.swap();
        aac.swap();
        acc.swap();
        ccc.swap();

        assert_eq!(pi.get_ref_pos_eager(&aaa), Some(aaa_pos));

        assert_eq!(pi.get_ref_pos_eager(&aac), Some(aac_pos));

        assert_eq!(pi.get_ref_pos_eager(&acc), Some(acc_pos));

        assert_eq!(pi.get_ref_pos_eager(&ccc), Some(ccc_pos));
    }

    #[test]
    fn validate_tiny() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(p).unwrap();
        pi.validate_self();
    }

    #[test]
    fn validate_tiny_sshash() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(p).unwrap();
        let unitig_set = pi.as_ref().clone();
        let sshash =
            SSHash::from_unitig_set_no_skew_index(unitig_set, 2, WyHashState::default()).unwrap();

        let pi = ModIndex::from_parts(
            pi.base.clone(),
            sshash,
            pi.as_u2pos().clone(),
            pi.as_refseqs().clone(),
        );
        pi.validate_self();
    }

    #[test]
    fn validate_tiny_refs_index() {
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(p).unwrap();
        pi.validate_self();
    }

    #[test]
    fn validate_small_txome_dense() {
        let p = to_abs_path(SMALL_TXOME_DENSE_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(p).unwrap();
        // pi.validate();
        pi.validate_self();
    }

    #[test]
    fn validate_yeast_dense() {
        let p = to_abs_path(YEAST_CHR01_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(p).unwrap();
        // pi.validate();
        pi.validate_self();
    }

    #[test]
    fn validate_yeast_sshash() {
        let p = to_abs_path(YEAST_CHR01_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(p).unwrap();
        let unitig_set = pi.as_ref().clone();
        let sshash = SSHash::from_unitig_set(unitig_set, 15, 32, WyHashState::default()).unwrap();

        let pi = ModIndex::from_parts(
            pi.base.clone(),
            sshash,
            pi.as_u2pos().clone(),
            pi.as_refseqs().clone(),
        );
        pi.validate_self();
    }
}
