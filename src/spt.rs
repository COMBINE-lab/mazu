use crate::cuttlefish::{CfFiles, CfSeqToken};
use crate::err::Result;
use crate::refseq::RefSeqCollection;
use crate::unitig_set::UnitigSet;
use crate::UnitigOcc;
use std::path::Path;

/// A spectrum preserving tiling with u64 packed unitig occurrence positions
///
/// Unitig occurrences are represented as an inverted list
// Note:
// - Perhaps update to also store polyN positions?
// - Could dump CuttlefishID-to-UnitigID mapping for future builds?
#[derive(Debug, Clone)]
pub struct SPT {
    pub unitigs: UnitigSet,
    pub ref_names: Vec<String>,
    pub ctable: Vec<u64>, // inv list of unitig occs (encoded into u64) // TODO rename this
    pub offsets: Vec<usize>, // prefix sum of offsets of inv lists
    pub ref_lens: Vec<usize>,
}

impl SPT {
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
        let word = self.ctable[i];
        UnitigOcc::decode_pf1(word)
    }
    pub fn get_occ(&self, uid: usize, i: usize) -> UnitigOcc {
        let ptr = self.offsets[uid] + i;
        let word = self.ctable[ptr];
        UnitigOcc::decode_pf1(word)
    }

    pub fn from_cf_prefix<P: AsRef<Path>>(fp_prefix: P) -> Result<Self> {
        let cf_files = CfFiles::new(fp_prefix);
        Self::from_cf_reduced_gfa(&cf_files)
    }

    pub fn from_cf_reduced_gfa(cf_files: &CfFiles) -> Result<Self> {
        let (unitigs, cfid2uid) = UnitigSet::from_cf_reduced_gfa(cf_files)?;
        let k = unitigs.k();
        let n_unitigs = unitigs.n_unitigs();

        let mut ufreq = vec![0; n_unitigs];
        let iter = cf_files.iter_tiling()?;

        // 1) Accumulate unitig occurence freqs to get ptrs into positions to write in inv list
        let mut ref_names: Vec<String> = Vec::new();
        for (ref_name, tiling) in iter {
            ref_names.push(ref_name);
            for token in tiling {
                if let CfSeqToken::Unitig { id, .. } = token {
                    let id = cfid2uid[&id];
                    ufreq[id] += 1
                }
            }
        }

        // 2) Prefix sum
        let offsets = crate::util::prefix_sum(&ufreq);
        let n_occs = offsets[n_unitigs];

        // 3) allocate inverted list
        let mut ctable_ptrs = offsets.clone(); // pointer to next word
        let mut ctable = vec![0; n_occs];
        let mut ref_lens = Vec::with_capacity(ref_names.len());

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
                        let word = occ.encode_pf1();

                        let ptr = ctable_ptrs[id]; // get ptr to position to insert
                        ctable[ptr] = word; // insert encoded word
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

        let spt = SPT::from_cf_prefix(TINY_CF_PREFIX).unwrap();

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
