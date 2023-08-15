use simple_sds::{int_vector::IntVector, ops::Pack};

use crate::{
    cuttlefish::CfFiles,
    index::{DenseUnitigTable, ModIndex, PiscemUnitigTable},
    kphf::{pfhash::PFHashDefault, sshash::SSHashDefault},
    spt::SPT,
    BaseIndex, ModIndexType, Result,
};

use std::path::Path;

pub type PufferfishDenseIndexDefault = ModIndex<PFHashDefault, DenseUnitigTable>;
pub type PiscemIndexDefault = ModIndex<SSHashDefault, PiscemUnitigTable>;

impl PufferfishDenseIndexDefault {
    pub fn from_cf_prefix<P: AsRef<Path>>(prefix: P) -> Result<Self> {
        let files = CfFiles::new(prefix);
        Self::from_cf_reduced_gfa(&files)
    }

    pub fn from_spt(spt: SPT) -> Result<Self> {
        let base = BaseIndex::new().set_index_type(ModIndexType::PufferfishDense);

        let refs = spt.get_ref_seq_collection();

        let sshash = PFHashDefault::from_unitig_set(spt.unitigs);
        let mut contig_offsets = IntVector::from(spt.offsets);
        contig_offsets.pack();

        let utab = DenseUnitigTable {
            ctable: spt.ctable,
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

    pub fn from_cf_reduced_gfa(cf_files: &CfFiles) -> Result<Self> {
        let spt = SPT::from_cf_reduced_gfa(cf_files)?;
        Self::from_spt(spt)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::index::Validate;
    use crate::test_utils::*;

    #[test]
    fn pufferfish_dense() {
        let prefix = YEAST_CF_PREFIX;
        let fasta = YEAST_CF_FA;

        let index = PufferfishDenseIndexDefault::from_cf_prefix(prefix).unwrap();
        index.validate_fasta(fasta);

        let mut index = index.as_streaming();
        index.validate_fasta(fasta);
    }
}
