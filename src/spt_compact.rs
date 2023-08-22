use simple_sds::int_vector::IntVector;
use simple_sds::ops::{Access, Vector};

use crate::cuttlefish::{CfFiles, CfSeqToken};
use crate::err::{Error, Result};
use crate::refseq::RefSeqCollection;
use crate::unitig_set::UnitigSet;
use crate::{Orientation, UnitigOcc};
// use std::io::{Error, ErrorKind};
use std::path::Path;

/// The mask we should apply to
/// obtain the position information
/// from a packed occurrence
/// (after right shifting 1 bit for the orientation)
static POS_MASKS: &[u64] = &[
    0x0,
    0x1,
    0x3,
    0x7,
    0xf,
    0x1f,
    0x3f,
    0x7f,
    0xff,
    0x1ff,
    0x3ff,
    0x7ff,
    0xfff,
    0x1fff,
    0x3fff,
    0x7fff,
    0xffff,
    0x1ffff,
    0x3ffff,
    0x7ffff,
    0xfffff,
    0x1fffff,
    0x3fffff,
    0x7fffff,
    0xffffff,
    0x1ffffff,
    0x3ffffff,
    0x7ffffff,
    0xfffffff,
    0x1fffffff,
    0x3fffffff,
    0x7fffffff,
    0xffffffff,
    0x1ffffffff,
    0x3ffffffff,
    0x7ffffffff,
    0xfffffffff,
    0x1fffffffff,
    0x3fffffffff,
    0x7fffffffff,
    0xffffffffff,
    0x1ffffffffff,
    0x3ffffffffff,
    0x7ffffffffff,
    0xfffffffffff,
    0x1fffffffffff,
    0x3fffffffffff,
    0x7fffffffffff,
    0xffffffffffff,
    0x1ffffffffffff,
    0x3ffffffffffff,
    0x7ffffffffffff,
    0xfffffffffffff,
    0x1fffffffffffff,
    0x3fffffffffffff,
    0x7fffffffffffff,
    0xffffffffffffff,
    0x1ffffffffffffff,
    0x3ffffffffffffff,
    0x7ffffffffffffff,
    0xfffffffffffffff,
    0x1fffffffffffffff,
    0x3fffffffffffffff,
    0x7fffffffffffffff,
];

impl UnitigOcc {
    /// Using the provided widths for the reference id
    /// and positions, encode this UnitigOcc as a u64
    /// NOTE: This assumes that the sum of the requried number
    /// of bits for the reference id, reference position, and orientation (1 bit)
    /// is <= 64.
    #[inline]
    pub fn encode_piscem(&self, ref_shift: usize) -> u64 {
        // debug_assert!(self.pos as u64 == (self.pos as u64 & pos_mask));
        let mut e = self.ref_id as u64;
        e <<= ref_shift;
        e |= (self.pos as u64) << 1;
        e |= if self.is_forward() { 1 } else { 0 };
        e
    }

    /// Given the provided shift for the reference id and the mask to extract
    /// the position, decode the provided `encoded_occ` into a UnitigOcc.
    pub fn decode_piscem(ref_shift: usize, pos_mask: u64, encoded_occ: u64) -> Self {
        UnitigOcc {
            ref_id: (encoded_occ >> ref_shift) as usize,
            pos: ((encoded_occ >> 1) & pos_mask) as usize,
            o: if encoded_occ & 0x1 == 1 {
                Orientation::Forward
            } else {
                Orientation::Backward
            },
        }
    }
}

#[derive(Debug, Clone)]
pub struct TileOccTable {
    pub(crate) ctable: IntVector, // inv list of unitig occs (encoded with the min. # of bits) // TODO rename this
    pub(crate) ref_shift: usize,  // the amount we have to shift an entry to get the reference
    pub(crate) pos_mask: u64, // the mask we apply to get the position after shifting off the orientation
}

impl TileOccTable {
    /// Create a new `TileOccTable` that will hold the relevant occurrence information.
    /// The input arguments allow the calculation of the number of bits required to hold
    /// each occ (roughly log2(max_ref_len) for reference position, log2(num_refs) for
    /// reference id, and 1 for the orientation). The `num_occs` is the total number of
    /// occurrences this `TileOccTable` will have to hold. Knowing this information is
    /// useful so that the memory can be allocated at once and there is no need to
    /// grow the `TileOccTable`.
    ///
    /// Returns `Ok(TileOccTable)` if the table could be created succesfully or
    /// `Error` if the table could not be constructed.
    fn from_occ_info(max_ref_len: usize, num_refs: usize, num_occs: usize) -> Result<Self> {
        if let Ok((pos_bits, _ref_bits, total_bits)) = required_num_bits(max_ref_len, num_refs) {
            let ctable = IntVector::with_len(num_occs, total_bits as usize, 0)
                .expect("Could not create IntVector to hold tile occurrences.");
            let ref_shift = (pos_bits as u64) + 1; // +1 for orientation bits
            let pos_mask = POS_MASKS[pos_bits as usize];
            Ok(Self {
                ctable,
                ref_shift: ref_shift as usize,
                pos_mask,
            })
        } else {
            Err(Error::Other(
                "Could not construct tile occurrence table.".to_string(),
            ))
        }
    }

    // #[inline]
    // fn set_encoded(&mut self, idx: usize, word: u64) {
    //     self.ctable.set(idx, word);
    // }

    /// Set the occurrence at index `idx` to the tile
    /// occurrence `occ`.
    #[inline]
    fn set_occ(&mut self, idx: usize, occ: UnitigOcc) {
        let e = occ.encode_piscem(self.ref_shift);
        self.ctable.set(idx, e);
    }

    /// Get the tile occurrence at index `idx` in this table.
    #[inline]
    fn get_occ(&self, idx: usize) -> UnitigOcc {
        let encoded_occ = self.ctable.get(idx);
        UnitigOcc::decode_piscem(self.ref_shift, self.pos_mask, encoded_occ)
    }

    /// The number of occurrences encoded in this table.
    #[inline]
    fn len(&self) -> usize {
        self.ctable.len()
    }

    /// Returns `true` if this tile table is empty (has no occurrences)
    /// and `false` otherwise.
    #[allow(unused)]
    #[inline]
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// A spectrum preserving tiling with packed unitig occurrence positions
///
/// Unitig occurrences are represented as an inverted list
// Note:
// - Perhaps update to also store polyN positions?
// - Could dump CuttlefishID-to-UnitigID mapping for future builds?
#[derive(Debug, Clone)]
pub struct SPTCompact {
    pub unitigs: UnitigSet,
    pub ref_names: Vec<String>,
    pub ctable: TileOccTable, // inv list of unitig occs (packed into minimal number of bits) // TODO rename this
    pub offsets: Vec<usize>,  // prefix sum of offsets of inv lists
    pub ref_lens: Vec<usize>,
}

/// The number of bits required to represent an integer with
/// a maximum value of `n`.
#[allow(unused)]
fn bit_width(n: usize) -> usize {
    if n == 0 {
        0
    } else {
        (crate::util::msb(n as u64) + 1) as usize
    }
}

/// Figures out the fixed size number of bits needed for each
/// tile occurrence entry given:
/// 1) the length of the longest observed reference sequence
/// 2) the number of references that need to be encoded
/// The return value is None if the result cannot be packed
/// into a single 64-bit integer. Otherwise, the return value
/// is a tuple (x, y, z) where
/// x : the number of bits required to encode the position on a reference
/// y : the number of bits required to encode the reference id
/// z : the total number of bits required to encode an entry (including orientation)
fn required_num_bits(longest_ref: usize, num_refs: usize) -> Result<(u8, u8, u8)> {
    let pos_bits = (longest_ref)
        .checked_ilog2()
        .ok_or(Error::other("Could not obtain log2 of longest reference."))?
        + 1;
    let ref_bits = (num_refs).checked_ilog2().ok_or(Error::other(
        "Could not obtain log2 of number of references.",
    ))? + 1;
    let ori_bits = 1_u32;
    let total_bits = ori_bits
        .checked_add(ref_bits)
        .ok_or(Error::other("Could not obtain required number of bits."))?
        .checked_add(pos_bits)
        .ok_or(Error::other("Could not obtain required number of bits."))?;
    if total_bits <= (u8::MAX as u32) {
        Ok((pos_bits as u8, ref_bits as u8, total_bits as u8))
    } else {
        Err(Error::other(
            "Could not compute number of bits required for each occ entry",
        ))
    }
}

impl SPTCompact {
    pub fn unitigs(&self) -> &UnitigSet {
        &self.unitigs
    }

    pub fn n_unitigs(&self) -> usize {
        self.unitigs.n_unitigs()
    }

    pub fn k(&self) -> usize {
        self.unitigs.k()
    }

    pub fn num_refs(&self) -> usize {
        self.ref_lens.len()
    }

    pub fn ref_name(&self, ri: usize) -> &str {
        &self.ref_names[ri]
    }

    pub fn ref_len(&self, ri: usize) -> usize {
        self.ref_lens[ri]
    }

    pub fn n_total_occs(&self) -> usize {
        self.ctable.len()
    }

    pub fn get_global_occ(&self, i: usize) -> UnitigOcc {
        self.ctable.get_occ(i)
    }

    pub fn get_occ(&self, uid: usize, i: usize) -> UnitigOcc {
        let ptr = self.offsets[uid] + i;
        self.ctable.get_occ(ptr)
    }

    pub fn from_cf_prefix<P: AsRef<Path>>(fp_prefix: P) -> Result<Self> {
        let cf_files = CfFiles::new(fp_prefix);
        Self::from_cf_reduced_gfa(&cf_files)
    }

    pub fn from_cf_reduced_gfa(cf_files: &CfFiles) -> Result<Self> {
        let (unitigs, cfid2uid) = UnitigSet::from_cf_reduced_gfa(cf_files)?;
        let k = unitigs.k();
        let n_unitigs = unitigs.n_unitigs();

        let mut ufreq: Vec<usize> = vec![0; n_unitigs];
        let iter = cf_files.iter_tiling()?;

        // 1) Accumulate unitig occurence freqs to get ptrs into positions to write in inv list
        // also record here, the maximum observed reference length
        let mut max_ref_len = 0_usize;
        let mut ref_names: Vec<String> = Vec::new();
        for (ref_name, tiling) in iter {
            ref_names.push(ref_name);
            let mut prev_was_unitig = false;
            let mut pos = 0; // start position of next tile
            for token in tiling {
                match token {
                    CfSeqToken::PolyN(n) => {
                        pos += n;
                        if prev_was_unitig {
                            pos += k - 1;
                        }
                        prev_was_unitig = false;
                    }

                    CfSeqToken::Unitig { id, .. } => {
                        let id = cfid2uid[&id];
                        let len = unitigs.unitig_len(id);
                        ufreq[id] += 1;
                        pos += len - k + 1;
                        prev_was_unitig = true;
                    }
                }
            }
            let len = if prev_was_unitig {
                pos + k - 1 // last tiling token was unitig
            } else {
                pos // last tile was polyN
            };
            max_ref_len = max_ref_len.max(len);
        }

        // 2) Prefix sum
        let offsets = crate::util::prefix_sum(&ufreq);
        let n_occs = offsets[n_unitigs];

        // 3) allocate inverted list
        let mut ctable_ptrs = offsets.clone(); // pointer to next word
        let mut ref_lens = Vec::with_capacity(ref_names.len());
        // compute the number of bits requried for each occurence
        // entry.
        // let bits_per_occ = required_num_bits(max_ref_len, ref_names.len())?;
        // TODO: @theJasonFan Right now the Error types for different Results are difficult to
        // satisfy. I propose we move to `anyhow` for handing / propagating results
        // or `thiserror` if we feel we need to specify different types for all errors
        // because they will propagate relevant information.
        let mut ctable = TileOccTable::from_occ_info(max_ref_len, ref_names.len(), n_occs)?;

        // 4) Iterate over tiling and insert occurrence info (encoded u64 words)
        let iter = cf_files.iter_tiling()?;
        for (ref_id, (_ref_name, tiling)) in iter.enumerate() {
            let mut prev_was_unitig = false;
            let mut pos = 0; // start position of next tile
            for token in tiling {
                match token {
                    CfSeqToken::PolyN(n) => {
                        pos += n; // skip polyN gap.
                        if prev_was_unitig {
                            pos += k - 1;
                        }
                        prev_was_unitig = false;
                    }
                    CfSeqToken::Unitig { id, o } => {
                        let id = cfid2uid[&id];
                        let len = unitigs.unitig_len(id);
                        let occ = UnitigOcc { ref_id, pos, o };

                        let ptr = ctable_ptrs[id]; // get ptr to position to insert
                        ctable.set_occ(ptr, occ); // insert encoded word
                        ctable_ptrs[id] += 1; // increment ptr by 1
                        pos += len - k + 1;
                        prev_was_unitig = true;
                    }
                }
            }

            let len = if prev_was_unitig {
                pos + k - 1 // last tiling token was unitig
            } else {
                pos // last tile was polyN
            };
            ref_lens.push(len); // push reference length.
        }

        Ok(Self {
            unitigs,
            ref_names,
            ctable,
            offsets,
            ref_lens,
        })
    }

    pub fn get_ref_seq_collection(&self) -> RefSeqCollection {
        RefSeqCollection {
            seq: None,
            prefix_sum: crate::util::prefix_sum(&self.ref_lens),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::test_utils::*;
    use crate::UnitigOcc;

    use crate::spt::SPT as SPT_;

    #[test]
    fn test_masks() {
        for i in 0..64 {
            let mask = POS_MASKS[i];
            assert_eq!(i, bit_width(mask as usize), "failed at i={}", i);
        }
    }

    #[test]
    fn compare_spt_impls() {
        let spt = SPTCompact::from_cf_prefix(TINY_CF_PREFIX).unwrap();
        let spt_old = SPT_::from_cf_prefix(TINY_CF_PREFIX).unwrap();
        assert_eq!(spt.offsets, spt_old.offsets);
        let n_occs = spt.n_total_occs();
        for i in 0..n_occs {
            let occ = spt.get_global_occ(i);
            let occ_ = spt_old.get_global_occ(i);
            assert_eq!(occ, occ_);
        }
        // panic!();

        let spt = SPTCompact::from_cf_prefix(YEAST_CF_PREFIX).unwrap();
        let spt_old = SPT_::from_cf_prefix(YEAST_CF_PREFIX).unwrap();
        assert_eq!(spt.offsets, spt_old.offsets);

        let n_occs = spt.n_total_occs();
        for i in 0..n_occs {
            let occ = spt.get_global_occ(i);
            let occ_ = spt_old.get_global_occ(i);
            assert_eq!(occ, occ_, "{}", i);
        }
    }

    #[test]
    fn tiny() {
        // tiny example with two seqs that have polyNs and are RC of eachother
        // tiny.seq:
        // Reference:1_Sequence:I	N3 3+ N1 2-
        // Reference:2_Sequence:I	2+ N1 3- N3
        //
        // tiny.seg:
        // 3	CACACACCAC
        // 2	CCTCAATACG
        //
        // Unitig IDs from cuttlefish is the minimum hash of kmers
        // which are re-ordered top to bottom in .seg file with ids 0-N.

        let spt = SPTCompact::from_cf_prefix(TINY_CF_PREFIX).unwrap();

        assert_eq!(spt.num_refs(), 2);
        assert_eq!(spt.n_total_occs(), 4);

        assert_eq!(spt.ref_name(0), "Reference:1_Sequence:I");
        assert_eq!(spt.ref_len(0), 3 + 10 + 1 + 10);

        // 0th occ of unitig 0
        let u0_occ0 = UnitigOcc {
            ref_id: 0,
            pos: 3,
            o: crate::Orientation::Forward,
        };
        assert_eq!(spt.get_occ(0, 0), u0_occ0);

        // 1th occ of unitig 0
        let u0_occ1 = UnitigOcc {
            ref_id: 1,
            pos: 10 + 1,
            o: crate::Orientation::Backward,
        };
        assert_eq!(spt.get_occ(0, 1), u0_occ1);

        // 0th occ of unitig 1
        let u1_occ0 = UnitigOcc {
            ref_id: 0,
            pos: 3 + 10 + 1,
            o: crate::Orientation::Backward,
        };

        assert_eq!(spt.get_occ(1, 0), u1_occ0);

        // 1th occ of unitig 1
        let u1_occ1 = UnitigOcc {
            ref_id: 1,
            pos: 0,
            o: crate::Orientation::Forward,
        };

        assert_eq!(spt.get_occ(1, 1), u1_occ1);
    }
}

#[cfg(test)]
mod test_spt_new {
    use crate::{Orientation, UnitigOcc};

    use super::POS_MASKS;

    #[test]
    fn occ_encoding() {
        let occ = UnitigOcc {
            ref_id: 0,
            pos: 1,
            o: Orientation::Backward,
        };

        let ref_shift = 3;
        let pos_mask = POS_MASKS[ref_shift - 1];

        let word = occ.encode_piscem(ref_shift);
        assert_eq!(word, 0b010);

        let occ_ = UnitigOcc::decode_piscem(ref_shift, pos_mask, word);
        assert_eq!(occ, occ_);
    }
}
