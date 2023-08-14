use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::{env, fmt, process};

// Yeast CHR01
#[allow(unused)]
pub const YEAST_CHR01_INDEX: &str = "test_data/pf1/yeast_chr01_index";

// Single contig
#[allow(unused)]
pub const TINY_INDEX: &str = "test_data/pf1/tiny_index";
#[allow(unused)]
pub const TINY_RC_INDEX: &str = "test_data/pf1/tiny-rc_index";

// Multi reference
#[allow(unused)]
pub const TINY_REFS_INDEX: &str = "test_data/pf1/tiny-multi-refs/tiny-multi-refs_index";
#[allow(unused)]
pub const TINY_REFS_FASTA: &str = "test_data/pf1/tiny-multi-refs/tiny-multi-refs.fasta";

// Small txome
#[allow(unused)]
pub const SMALL_TXOME_DENSE_INDEX: &str = "test_data/pf1/small_txome_index";
#[allow(unused)]
pub const SMALL_TXOME_SPARSE_INDEX: &str = "test_data/pf1/small_txome_index_sparse";

#[allow(unused)]
pub fn assert_path_eq(p: PathBuf, s: &str) {
    assert_eq!(p.to_str().unwrap(), s);
}

#[allow(unused)]
pub fn to_abs_path<P: AsRef<Path>>(path: P) -> PathBuf {
    let workdir = env!("CARGO_MANIFEST_DIR");
    let workdir = Path::new(&workdir);
    workdir.join(path)
}

/// Temp file generateion taken from simple_sds
// Counter used for temporary file names.
#[allow(unused)]
static TEMP_FILE_COUNTER: AtomicUsize = AtomicUsize::new(0);

/// Returns a name for a temporary file using the provided name part.
///
/// # Examples
///
/// ```
/// use simple_sds::serialize;
///
/// let filename = serialize::temp_file_name("example");
/// assert!(filename.into_os_string().into_string().unwrap().contains("example"));
/// ```
#[allow(unused)]
pub fn temp_file_name(name_part: &str) -> PathBuf {
    let count = TEMP_FILE_COUNTER.fetch_add(1, Ordering::SeqCst);
    let mut buf = env::temp_dir();
    buf.push(format!("{}_{}_{}", name_part, process::id(), count));
    buf
}

/// Taken from itertools crate
/// https://docs.rs/itertools/0.5.9/src/itertools/lib.rs.html#1576-1599
///
/// Assert that two iterators produce equal sequences, with the same
/// semantics as *equal(a, b)*.
///
/// **Panics** on assertion failure with a message that shows the
/// two iteration elements.
///
/// ```ignore
/// assert_equal("exceed".split('c'), "excess".split('c'));
/// // ^PANIC: panicked at 'Failed assertion Some("eed") == Some("ess") for iteration 1',
/// ```
#[allow(unused)]
pub fn assert_iter_eq<I, J>(a: I, b: J)
where
    I: IntoIterator,
    J: IntoIterator,
    I::Item: fmt::Debug + PartialEq<J::Item>,
    J::Item: fmt::Debug,
{
    let mut ia = a.into_iter();
    let mut ib = b.into_iter();
    let mut i = 0;
    loop {
        match (ia.next(), ib.next()) {
            (None, None) => return,
            (Some(a), Some(b)) => {
                assert_eq!(a, b, "\nFailed for iteration {i}\n",);
                i += 1;
            }
            (a, b) => {
                panic!("\nFailed for iteration {i}\n{a:#?}\n=/=\n{b:#?}\n");
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    #[should_panic]
    fn not_iter_eq() {
        assert_iter_eq([0, 1], [0, 0]);
    }

    #[test]
    fn iter_eq() {
        assert_iter_eq([0, 1], [0, 1]);
    }
}
