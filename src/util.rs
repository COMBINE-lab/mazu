use std::{
    cell::UnsafeCell,
    io::{BufRead, Lines},
    iter::Peekable,
    ops::Add,
};

/// Convenience trait to get the "zero" (additive identity) of a numeric type
pub trait Zero {
    fn zero() -> Self;
}

impl Zero for u64 {
    fn zero() -> Self {
        0
    }
}

impl Zero for usize {
    fn zero() -> Self {
        0
    }
}

/// Compute the prefix sum of a slice, returning a vector
pub fn prefix_sum<T: Add<Output = T> + Zero + Copy>(xs: &[T]) -> Vec<T>
where
{
    let mut accum = T::zero();
    let mut res = Vec::with_capacity(xs.len() + 1);
    for &x in xs {
        res.push(accum);
        accum = accum + x;
    }
    res.push(accum);

    res
}

/// Return index of most significant bit
///
/// ```
/// use mazu::util::msb;
/// assert_eq!(msb(0), 0);
/// assert_eq!(msb(0b01), 0);
/// assert_eq!(msb(0b0101), 2);
/// ```
pub fn msb(n: u64) -> u64 {
    // https://github.com/tomarrell/rust-elias-fano/blob/master/src/utils.rs
    if n == 0 {
        0
    } else {
        ((u64::BITS as u64) - 1) - n.leading_zeros() as u64
    }
}

// Unsafe slice allowing completely unsafe writes to `slice[i]` via raw pointers.
// Useful for parallel construction. Use with caution.
// https://stackoverflow.com/questions/55939552/simultaneous-mutable-access-to-arbitrary-indices-of-a-large-vector-that-are-guar
#[derive(Copy, Clone)]
pub(crate) struct UnsafeSlice<'a, T> {
    slice: &'a [UnsafeCell<T>],
}
unsafe impl<'a, T: Send + Sync> Send for UnsafeSlice<'a, T> {}
unsafe impl<'a, T: Send + Sync> Sync for UnsafeSlice<'a, T> {}

impl<'a, T> UnsafeSlice<'a, T> {
    pub fn new(slice: &'a mut [T]) -> Self {
        let ptr = slice as *mut [T] as *const [UnsafeCell<T>];
        Self {
            slice: unsafe { &*ptr },
        }
    }

    /// # Safety
    ///
    /// It is UB if two threads write to the same index without
    /// synchronization. This function itself does not protect the
    /// underlying resource, so it is the responsibility of the
    /// caller to ensure that multiple threads do not attempt to
    /// write to the same index simultaneously.
    pub unsafe fn write(&self, i: usize, value: T) {
        let ptr = self.slice[i].get();
        *ptr = value;
    }
}

/// Simple reader for a FASTA file. Implements [Iterator] over [FastaRecord]s.
///
/// Returned records does no processing to FASTA record names, returning
/// FASTA header verbatim (following `>`).
pub struct FastaReader<T: BufRead> {
    lines: Peekable<Lines<T>>,
    seq_buf: String,
}

// A FASTA record.
#[derive(Debug, PartialEq, Eq)]
pub struct FastaRecord {
    pub name: String,
    pub seq: String,
}

impl AsRef<[u8]> for FastaRecord {
    fn as_ref(&self) -> &[u8] {
        self.seq.as_bytes()
    }
}

impl<T: BufRead> FastaReader<T> {
    pub fn new(buf_read: T) -> Self {
        let lines = buf_read.lines().peekable();
        Self {
            lines,
            seq_buf: String::new(),
        }
    }

    fn should_flush(&mut self) -> bool {
        let line = self.lines.peek();
        if let Some(rline) = line {
            let line = rline.as_ref().expect("could not parse line");
            line.starts_with('>')
        } else {
            true
        }
    }
}

impl<T: BufRead> Iterator for FastaReader<T> {
    type Item = FastaRecord;

    fn next(&mut self) -> Option<Self::Item> {
        let name = self.lines.next()?.expect("could not parse line");
        let (_, name) = name.split_at(1);
        self.seq_buf.clear();

        while !self.should_flush() {
            let seq = self.lines.next().unwrap().unwrap();
            self.seq_buf.push_str(&seq);
        }

        Some(FastaRecord {
            name: name.to_string(),
            seq: self.seq_buf.clone(),
        })
    }
}

#[cfg(test)]
mod test {
    use std::io::BufReader;

    use super::*;

    #[test]
    fn fasta_reader() {
        let input = ">A\n\
        ACG\n\
        AC\n\
        >B\n\
        ACG\n\
        ACG\n\
        C"
        .as_bytes();
        let rdr = BufReader::new(input);
        let rdr = FastaReader::new(rdr);

        let v = Vec::from_iter(rdr);
        assert_eq!(
            v[0],
            FastaRecord {
                name: "A".to_string(),
                seq: "ACGAC".to_string()
            }
        );
        assert_eq!(
            v[1],
            FastaRecord {
                name: "B".to_string(),
                seq: "ACGACGC".to_string()
            }
        );
    }
}
