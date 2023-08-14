pub mod cuttlefish;
pub mod elias_fano;
mod err;
pub mod kphf;
pub mod pf1;
pub mod refseq;
pub mod spt;
pub mod spt_compact;
pub mod unitig_set;
pub mod util;
pub mod wm;

#[cfg(test)]
pub mod test_utils;

// pub mod cc_index;
pub mod index;
pub use index::*;

pub use err::{Error, Result};

#[allow(non_camel_case_types)]
type km_size_t = usize;

/// Orientation `enum` for sequences (unitigs, k-mers, etc.,)
#[derive(Copy, Debug, Clone, PartialEq, Eq)]
pub enum Orientation {
    Forward,
    Backward,
}

impl Orientation {
    pub fn reverse(self) -> Self {
        match self {
            Self::Forward => Self::Backward,
            Self::Backward => Self::Forward,
        }
    }
}

impl Orientation {
    /// Sensible conversions to/from unsigned integers
    pub fn to_u8(self) -> u8 {
        match self {
            Orientation::Forward => 1,
            Orientation::Backward => 0,
        }
    }

    /// Sensible conversions to/from unsigned integers
    pub fn to_u64(self) -> u64 {
        match self {
            Orientation::Forward => 1,
            Orientation::Backward => 0,
        }
    }

    /// Sensible conversions to/from unsigned integers
    ///
    /// Returns `Orientation::Forward` if `w > 0`
    pub fn from_u8(w: u8) -> Self {
        if w > 0 {
            Orientation::Forward
        } else {
            Orientation::Backward
        }
    }

    /// Sensible conversions to/from unsigned integers
    pub fn from_u64(w: u64) -> Self {
        if w > 0 {
            Orientation::Forward
        } else {
            Orientation::Backward
        }
    }
}
