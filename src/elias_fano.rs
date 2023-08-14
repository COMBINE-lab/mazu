use simple_sds::{
    bit_vector::BitVector,
    int_vector::IntVector,
    ops::{Access, Push, Select, Vector},
    raw_vector::{AccessRaw, RawVector},
};

use crate::err::{Error, Result};

use serde::{Deserialize, Serialize};

/// Elias-Fano compressed representation for monotonic sequence
///
/// ```
/// use mazu::elias_fano::EFVector;
/// let xs = [5, 8, 8, 15, 32]; // Example from Vigna 2013
/// let ys = EFVector::from_iter(xs.clone().into_iter(), 32).unwrap();
/// for (i, &x) in xs.iter().enumerate() {
///     assert_eq!(ys.get(i), x);
/// }
/// ```
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct EFVector {
    u: u64, // universe size
    // n: u64, // length
    l: u64, // low_bit_width
    high_bits: BitVector,
    low_bits: IntVector,
}

impl EFVector {
    pub fn num_bits(&self) -> usize {
        std::mem::size_of::<u64>() * 3 * 8 + self.high_bits.num_bits() + self.low_bits.num_bits()
    }

    pub fn len(&self) -> usize {
        self.low_bits.len()
    }

    pub fn is_empty(&self) -> bool {
        self.low_bits.len() == 0
    }

    pub fn u(&self) -> usize {
        self.u as usize
    }

    pub fn low_bit_width(&self) -> usize {
        self.l as usize
    }

    pub fn from_iter<I>(it: I, u: usize) -> crate::err::Result<Self>
    where
        I: Iterator<Item = u64> + ExactSizeIterator,
    {
        assert!(u <= (u64::MAX as usize));
        let n = it.len() as u64;
        let u = u as u64;
        let l = crate::util::msb(u / n);

        // hack to avoid l == 0
        let l = if l == 0 {
            log::warn!(
                "Elias-Fano lower-bit-width computed to be 0 (from u={}, n={}), setting l=1 instead",
                u,
                n
            );
            1
        } else {
            l
        };

        let low_mask = (1 << l) - 1;

        let mut low_bits = IntVector::with_capacity(it.len(), l as usize).unwrap();

        let hb_len = n + (u >> l); // total gap is u >> l, n stop bits,
        let mut high_bits = RawVector::with_len(hb_len as usize, false);

        let mut last = 0;
        let mut h0s = 0; // # of 0s preceeding high bits for elem i (stop bits)

        for (i, x) in it.enumerate() {
            if x < last {
                return Err(Error::EFNotMonotone);
            }

            let lbits = x & low_mask;
            low_bits.push(lbits);

            let hbit_gap = ((x >> l) - (last >> l)) as usize;

            // offset is # of preceeding 0s + 1 bits + gap.
            let set_bit = h0s + i + hbit_gap;
            high_bits.set_bit(set_bit, true);

            last = x;
            h0s += hbit_gap;
        }

        let mut high_bits = BitVector::from(high_bits);
        high_bits.enable_select();

        Ok(Self {
            //n,
            u,
            l,
            low_bits,
            high_bits,
        })
    }

    pub fn get(&self, i: usize) -> u64 {
        let high = self.high_bits.select(i).unwrap() - i;
        let high = high as u64;
        let low = self.low_bits.get(i);

        (high << self.l) | low
    }

    pub fn from_slice(xs: &[u64]) -> Result<Self> {
        if xs.is_empty() {
            Err(Error::EFEmpty)
        } else {
            let u = xs[xs.len() - 1];
            Self::from_iter(xs.iter().copied(), u as usize)
        }
    }

    pub fn from_usize_slice(xs: &[usize]) -> Result<Self> {
        if xs.is_empty() {
            Err(Error::EFEmpty)
        } else {
            let u = xs[xs.len() - 1];
            Self::from_iter(xs.iter().map(|x| *x as u64), u)
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn vigna_fig1() {
        let xs = [5, 8, 8, 15, 32];

        let ys = EFVector::from_iter(xs.clone().into_iter(), 32).unwrap();

        for (i, &x) in xs.iter().enumerate() {
            assert_eq!(ys.get(i), x);
        }
    }

    #[should_panic]
    #[test]
    fn not_monotone() {
        let xs = [5, 8, 7, 15, 32];

        EFVector::from_iter(xs.clone().into_iter(), 32).expect("Monotone seq");
    }
}
