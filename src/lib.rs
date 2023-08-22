// #[warn(missing_docs)]
#![warn(missing_debug_implementations)]

// Basic data structures, utilities, errors
/// Elias-Fano compression
pub mod elias_fano;
pub mod err;
pub mod util;
pub mod wm;

// Interaction with COMBINE-lab projects
pub mod cuttlefish;
pub mod pf1;

// Index building blocks
pub mod index;
pub mod kphf;
pub mod refseq;
pub mod spt;
pub mod spt_compact;
pub mod unitig_set;

#[cfg(test)]
pub mod test_utils;

// pub mod cc_index;
pub use index::*;

pub use err::{Error, Result};

pub fn get_mazu_version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

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

    /// Sensible conversions to/from unsigned integers
    pub fn to_u8(&self) -> u8 {
        match self {
            Orientation::Forward => 1,
            Orientation::Backward => 0,
        }
    }

    /// Sensible conversions to/from unsigned integers
    pub fn to_u64(&self) -> u64 {
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
