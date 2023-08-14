// All things to do with Cpp impl of pufferfish
use std::fs::File;
// use std::io::{Error, ErrorKind};
use std::io::Read;

use std::path::Path;

use simple_sds::{bit_vector::BitVector, int_vector::IntVector, raw_vector::RawVector};

use kmers::naive_impl::seq_vector::SeqVector;

use crate::{Error, Result};

/******************************************************************************/
// ReadFrom --
//   convenience trait to read bytes into types directly from file
//   pointers
/******************************************************************************/
pub trait ReadFrom {
    fn read_from(f: &mut dyn Read) -> Result<Self>
    where
        Self: Sized;
}

impl ReadFrom for u64 {
    fn read_from(f: &mut dyn Read) -> Result<Self> {
        let mut val: [u8; 8] = [0; 8];
        f.read_exact(&mut val)?;
        let val = u64::from_ne_bytes(val);
        Ok(val)
    }
}
impl ReadFrom for f64 {
    fn read_from(f: &mut dyn Read) -> Result<Self> {
        let mut val: [u8; 8] = [0; 8];
        f.read_exact(&mut val)?;
        let val = f64::from_ne_bytes(val);
        Ok(val)
    }
}

impl ReadFrom for u32 {
    fn read_from(f: &mut dyn Read) -> Result<Self> {
        let mut val: [u8; 4] = [0; 4];
        f.read_exact(&mut val)?;
        let val = u32::from_ne_bytes(val);
        Ok(val)
    }
}

impl ReadFrom for i32 {
    fn read_from(f: &mut dyn Read) -> Result<Self> {
        let mut val: [u8; 4] = [0; 4];
        f.read_exact(&mut val)?;
        let val = i32::from_ne_bytes(val);
        Ok(val)
    }
}

impl ReadFrom for Vec<u64> {
    fn read_from(f: &mut dyn Read) -> Result<Self> {
        type T = u64;
        let mut buf = Vec::new();
        let n_bytes = f.read_to_end(&mut buf)?;
        let sizeof_t = std::mem::size_of::<T>();

        if (n_bytes % sizeof_t) != 0 {
            let msg = format!("Cannot read {n_bytes} bytes not divisible by {sizeof_t}");
            return Err(Error::InvalidData(msg));
        }
        let n_words = n_bytes / sizeof_t;

        let mut v = Vec::with_capacity(n_words);

        for chunk in buf.chunks_exact(sizeof_t) {
            // unwrap is fine here. 8 bytes guaranteed by chunks_exact
            let chunk = chunk.try_into().unwrap();
            let word = T::from_ne_bytes(chunk);
            v.push(word);
        }
        Ok(v)
    }
}

impl ReadFrom for Vec<u32> {
    fn read_from(f: &mut dyn Read) -> Result<Self> {
        type T = u32;
        let mut buf = Vec::new();
        let n_bytes = f.read_to_end(&mut buf)?;
        let sizeof_t = std::mem::size_of::<T>();

        if (n_bytes % sizeof_t) != 0 {
            let msg = format!("Cannot read {n_bytes} bytes not divisible by {sizeof_t}");
            return Err(Error::InvalidData(msg));
        }
        let n_words = n_bytes / sizeof_t;

        let mut v = Vec::with_capacity(n_words);

        for chunk in buf.chunks_exact(sizeof_t) {
            // unwrap is fine here. 8 bytes guaranteed by chunks_exact
            let chunk = chunk.try_into().unwrap();
            let word = T::from_ne_bytes(chunk);
            v.push(word);
        }
        Ok(v)
    }
}

pub fn read_u64_vec_with_len(f: &mut dyn Read, len: usize) -> Result<Vec<u64>> {
    type T = u64;
    let n_bytes = len * std::mem::size_of::<T>();
    let mut f = f.take(n_bytes as u64);
    Vec::read_from(&mut f)
}

pub fn read_u32_vec_with_len(f: &mut dyn Read, len: usize) -> Result<Vec<u32>> {
    type T = u32;
    let n_bytes = len * std::mem::size_of::<T>();
    let mut f = f.take(n_bytes as u64);
    Vec::read_from(&mut f)
}

pub trait FromCereal {
    fn read_from_cereal_archive(f: &mut dyn Read) -> Result<Self>
    where
        Self: Sized;

    fn load_from_cereal_archive<P: AsRef<Path>>(p: P) -> Result<Self>
    where
        Self: Sized,
    {
        let mut f = File::open(p)?;
        Self::read_from_cereal_archive(&mut f)
    }
}

impl FromCereal for Vec<String> {
    // Deserialize Cereal input archive for std::vec of strings
    // <len_of_vec><n_chars><str>...<nchars><str>
    fn read_from_cereal_archive(f: &mut dyn Read) -> Result<Self> {
        let len = u64::read_from(f)? as usize;
        let mut v = Vec::with_capacity(len);
        for _ in 0..len {
            let n_chars = u64::read_from(f)? as usize;
            let mut buf = vec![0u8; n_chars];
            f.read_exact(&mut buf)?;

            let s = String::from_utf8(buf).unwrap();
            v.push(s);
        }
        Ok(v)
    }
}

impl FromCereal for Vec<u32> {
    // Deserialize Cereal input archive for std::vec
    // <len_of_vec><word>...<word>
    fn read_from_cereal_archive(f: &mut dyn Read) -> Result<Self> {
        let len = u64::read_from(f)? as usize;
        read_u32_vec_with_len(f, len)
    }
}

impl FromCereal for Vec<u64> {
    // Deserialize Cereal input archive for std::vec
    // <len_of_vec><word>...<word>
    fn read_from_cereal_archive(f: &mut dyn Read) -> Result<Self> {
        let len = u64::read_from(f)? as usize;
        read_u64_vec_with_len(f, len)
    }
}

pub trait DeserializeFromCpp
where
    Self: Sized,
{
    fn deserialize_from_cpp<P: AsRef<Path>>(p: P) -> Result<Self>;
}

#[inline]
#[allow(dead_code)]
pub fn get_bits_per_element<P: AsRef<Path>>(p: P) -> usize {
    let mut f = File::open(p).unwrap();
    let _static_flag = u64::read_from(&mut f).unwrap();
    let width = u64::read_from(&mut f).unwrap();
    width as usize
}

pub trait FromCompact {
    fn from_compact_serialized<P: AsRef<Path>>(p: P) -> Result<Self>
    where
        Self: Sized;
}

impl FromCompact for SeqVector {
    fn from_compact_serialized<P: AsRef<Path>>(p: P) -> Result<Self> {
        let iv = IntVector::from_compact_serialized(p)?;
        Ok(Self::from(iv))
    }
}

impl FromCompact for RawVector {
    fn from_compact_serialized<P: AsRef<Path>>(p: P) -> Result<Self> {
        let iv = IntVector::from_compact_serialized(p)?;
        Ok(RawVector::from(iv))
    }
}

impl FromCompact for BitVector {
    fn from_compact_serialized<P: AsRef<Path>>(p: P) -> Result<Self> {
        let rv = RawVector::from_compact_serialized(p)?;
        Ok(BitVector::from(rv))
    }
}

impl FromCompact for IntVector {
    fn from_compact_serialized<P: AsRef<Path>>(p: P) -> Result<Self> {
        let mut f = File::open(p)?;
        let _static_flag = u64::read_from(&mut f)?;
        // What is this? this just sets the bits per element
        //assert_eq!(static_flag, 0); // TODO

        let width = u64::read_from(&mut f)? as usize;
        assert!(width > 0);

        let len = u64::read_from(&mut f)? as usize;
        let _capacity = u64::read_from(&mut f)?; // we don't use this

        let words = Vec::read_from(&mut f)?;
        let rv = RawVector::from_parts(len * width, words);

        let iv = IntVector::from_parts(len, width, rv);

        Ok(iv)
    }
}

#[cfg(test)]
mod int_vec_test {
    use super::*;
    use simple_sds::ops::{Access, Vector};

    #[test]
    fn test() {
        let pos_fp = "test_data/pf1/tiny_index/pos.bin";
        let workdir = env!("CARGO_MANIFEST_DIR");
        let workdir = Path::new(workdir);
        let pos_fp = workdir.join(pos_fp);
        let iv = IntVector::from_compact_serialized(&pos_fp).unwrap();
        let width = get_bits_per_element(&pos_fp);
        assert_eq!(iv.width(), width);
        assert_eq!(iv.width(), 3);
        assert_eq!(iv.len(), 4);

        let correct_ints = vec![3, 2, 0, 1];
        for (loaded, correct) in iv.iter().zip(correct_ints.iter()) {
            assert_eq!(loaded, *correct);
        }
    }
}
#[cfg(test)]
mod read_width_tests {
    use super::*;

    #[test]
    fn test_get_bits_per_element_yeast_chr01() {
        let pos_fp = "test_data/pf1/yeast_chr01_index/pos.bin";
        // 221918 kmers in yeastchr01, 221918.log2.ceil()
        test_get_bits_per_element(pos_fp, 18);
    }

    #[test]
    fn test_get_bits_per_element_tiny() {
        let pos_fp = "test_data/pf1/tiny_index/pos.bin";
        // Note in line 614 of PufferfishIndexer.cpp:
        //     size_t w = std::log2(tlen) + 1;
        // so 4 kmers in tiny example means, width = 2 + 1
        test_get_bits_per_element(pos_fp, 3);
    }

    fn test_get_bits_per_element(pos_fp: &str, correct_w: usize) {
        let workdir = env!("CARGO_MANIFEST_DIR");
        let workdir = Path::new(workdir);
        let pos_fp = workdir.join(pos_fp);
        let loaded_w = get_bits_per_element(pos_fp);

        assert_eq!(loaded_w, correct_w);
    }
}
