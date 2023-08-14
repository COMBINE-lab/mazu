use simple_sds::int_vector::IntVector;
use std::{fs::File, io::Read, path::Path};

use super::cpp::{DeserializeFromCpp, FromCereal, FromCompact};
use super::PF1FilePaths;

use crate::index::DenseUnitigTable;
use crate::Result;

impl DeserializeFromCpp for DenseUnitigTable {
    fn deserialize_from_cpp<P: AsRef<Path>>(dpath: P) -> Result<Self> {
        let files = PF1FilePaths::new(dpath);
        let (ref_names, ref_exts, ctable) = Self::deserialize_from_cpp_helper(files.ctable)?;
        let contig_offsets = IntVector::from_compact_serialized(files.ctg_offsets)?;
        Ok(Self {
            ctable,
            contig_offsets,
            ref_names,
            _ref_exts: ref_exts,
        })
    }
}

impl DenseUnitigTable {
    /**************************************************************************/
    // Helpers
    /**************************************************************************/
    fn deserialize_from_cpp_helper<P: AsRef<Path>>(
        p: P,
    ) -> Result<(Vec<String>, Vec<u32>, Vec<u64>)> {
        // Helper fn to deserialize serialized bits from cpp impl
        // that stores info as a tuple of vectors
        let mut f = File::open(p)?;

        // Deserialize Cereal input archive for Vec of std::string
        let ref_names = Vec::read_from_cereal_archive(&mut f)?;

        // Deserialize Cereal input archive for Vec of uint32
        let ref_exts = Vec::read_from_cereal_archive(&mut f)?;

        // Deserialize Cereal input archive for Vec of uint64
        let ctable = Vec::read_from_cereal_archive(&mut f)?;

        let mut buf = Vec::new();
        let n_bytes = f.read_to_end(&mut buf).unwrap();
        assert_eq!(n_bytes, 0);

        Ok((ref_names, ref_exts, ctable))
    }
}
