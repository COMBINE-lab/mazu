use crate::{cuttlefish::CfFiles, unitig_set::UnitigSet};
use std::path::{Path, PathBuf};

pub const TINY_CF_PREFIX: &str = "test_data/cf/tiny/tiny";
pub const TINY_CF_FA: &str = "test_data/cf/tiny/tiny.fa";

// Note: with w=15, some minimizers occur > 32 times
// setting sshash skew param too high (e.g. 64) will result in no skew minimizers
pub const YEAST_CF_PREFIX: &str = "test_data/cf/yeast_chr7/yeast_chr7";
pub const YEAST_CF_FA: &str = "test_data/cf/yeast_chr7/yeast_chr7.fa";

pub fn to_abs_path<P: AsRef<Path>>(path: P) -> PathBuf {
    let workdir = env!("CARGO_MANIFEST_DIR");
    let workdir = Path::new(&workdir);
    workdir.join(path)
}

// TODO: this really should just be part of UnitigSet.
pub fn load_unitigs<P: AsRef<Path>>(prefix: P) -> UnitigSet {
    let files = CfFiles::new(prefix);
    let (unitigs, _) = UnitigSet::from_cf_reduced_gfa(&files).unwrap();
    unitigs
}
