use simple_sds::{int_vector::IntVector, ops::Pack};
use std::path::Path;

use crate::{
    cuttlefish::CfFiles,
    dense_unitig_table::PiscemUnitigTable,
    kphf::{sshash::SSHashDefault, WyHashState},
    BaseIndex, ModIndex, ModIndexType, Result,
};

use crate::spt_compact::SPTCompact;

type PiscemIndex = ModIndex<SSHashDefault, PiscemUnitigTable>;

impl PiscemIndex {
    pub fn as_sshash(&self) -> &SSHashDefault {
        &self.k2u
    }
    pub fn from_cf_prefix<P: AsRef<Path>>(prefix: P, w: usize, skew_param: usize) -> Result<Self> {
        let files = CfFiles::new(prefix);
        Self::from_cf_reduced_gfa(&files, w, skew_param)
    }

    pub fn from_spt(spt: SPTCompact, w: usize, skew_param: usize) -> Result<Self> {
        let base = BaseIndex::new().set_index_type(ModIndexType::Piscem);
        let refs = spt.get_ref_seq_collection();

        let bh = WyHashState::default();
        let sshash = SSHashDefault::from_unitig_set(spt.unitigs, w, skew_param, bh)?;
        let mut contig_offsets = IntVector::from(spt.offsets);
        contig_offsets.pack();

        let utab = PiscemUnitigTable {
            ctable: spt.ctable.ctable,
            ref_shift: spt.ctable.ref_shift,
            pos_mask: spt.ctable.pos_mask,
            contig_offsets,
            ref_names: spt.ref_names,
            _ref_exts: Vec::new(),
        };

        let index = ModIndex {
            base,
            u2pos: utab,
            k2u: sshash,
            refs,
        };

        Ok(index)
    }

    pub fn from_cf_reduced_gfa(cf_files: &CfFiles, w: usize, skew_param: usize) -> Result<Self> {
        let spt = SPTCompact::from_cf_reduced_gfa(cf_files)?;
        Self::from_spt(spt, w, skew_param)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{index::Validate, test_utils::*};

    #[test]
    fn validate_tiny() {
        let w = 3;
        let skew_param = 2;
        let prefix = TINY_CF_PREFIX;
        let fasta = TINY_CF_FA;

        let index = PiscemIndex::from_cf_prefix(prefix, w, skew_param).unwrap();
        index.validate_fasta(fasta);
    }

    #[test]
    fn validate_yeast() {
        let w = 15;
        let skew_param = 32;
        let prefix = YEAST_CF_PREFIX;
        let fasta = YEAST_CF_FA;

        let index = PiscemIndex::from_cf_prefix(prefix, w, skew_param).unwrap();
        index.validate_fasta(fasta);

        // just some sanity checks
        assert!(index.as_sshash().n_minimizers() < index.n_kmers());
        assert_ne!(index.as_sshash().n_kmers_in_skew_index(), 0); // make sure the test case is not trivial
    }

    #[test]
    fn validate_yeast_spt_new_streaming() {
        let w = 15;
        let skew_param = 32;
        let prefix = YEAST_CF_PREFIX;
        let fasta = YEAST_CF_FA;

        let index = PiscemIndex::from_cf_prefix(prefix, w, skew_param).unwrap();
        let mut index = index.as_streaming();
        index.validate_fasta(fasta);
    }
}
