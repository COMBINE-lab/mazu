use std::ops::Range;

use simple_sds::{
    int_vector::IntVector,
    ops::{Access, Vector},
};

use super::{U2Pos, UnitigOcc};

use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DenseUnitigTable {
    pub(crate) ctable: Vec<u64>,
    pub(crate) contig_offsets: IntVector,
    pub(crate) ref_names: Vec<String>, // TODO move this into refseq class
    pub(crate) _ref_exts: Vec<u32>,
}

impl DenseUnitigTable {
    pub fn get_encoded_contig_occs(&self, ctg_id: usize) -> &[u64] {
        let start_idx = self.contig_offsets.get(ctg_id) as usize;
        let end_idx = self.contig_offsets.get(ctg_id + 1) as usize;
        let r = start_idx..end_idx;
        return self.ctable.get(r).unwrap();
    }

    pub fn num_ctg_occs(&self, ctg_id: usize) -> usize {
        let s = self.contig_offsets.get(ctg_id) as usize;
        let e = self.contig_offsets.get(ctg_id + 1) as usize;
        e - s
    }

    pub fn num_total_occs(&self) -> usize {
        self.ctable.len()
    }

    pub fn num_contigs(&self) -> usize {
        self.contig_offsets.len() - 1
    }

    pub fn get_offsets(&self) -> &IntVector {
        &self.contig_offsets
    }

    pub fn get_ref_name(&self, id: usize) -> Option<&String> {
        self.ref_names.get(id)
    }

    pub fn get_ref_names(&self) -> &[String] {
        &self.ref_names
    }
}

impl U2Pos for DenseUnitigTable {
    type EncodedOccs<'a> = &'a [u64] where Self: 'a;

    fn encoded_unitig_occs(&self, ui: usize) -> Self::EncodedOccs<'_> {
        let start_idx = self.contig_offsets.get(ui) as usize;
        let end_idx = self.contig_offsets.get(ui + 1) as usize;
        let r = start_idx..end_idx;
        &self.ctable[r]
    }

    fn decode_unitig_occs<T>(&self, _k2u: &T, occs: Self::EncodedOccs<'_>) -> Vec<UnitigOcc> {
        occs.iter().map(|&enc| UnitigOcc::decode_pf1(enc)).collect()
    }

    fn n_unitig_occs(&self, occs: Self::EncodedOccs<'_>) -> usize {
        occs.len()
    }

    fn empty_encoded_occs(&self) -> Self::EncodedOccs<'_> {
        &[]
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OccsRange {
    start: usize,
    end: usize,
}

impl OccsRange {
    pub fn empty() -> Self {
        Self {
            start: usize::MAX,
            end: usize::MAX,
        }
    }

    pub fn as_range(&self) -> Range<usize> {
        self.start..self.end
    }

    pub fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }

    pub fn len(&self) -> usize {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PiscemUnitigTable {
    pub(crate) ref_shift: usize,
    pub(crate) pos_mask: u64,

    pub(crate) ctable: IntVector,
    pub(crate) contig_offsets: IntVector,
    pub(crate) ref_names: Vec<String>, // TODO move this into refseq class
    pub(crate) _ref_exts: Vec<u32>,
}

impl PiscemUnitigTable {
    #[inline]
    fn decode_occ_word(&self, word: u64) -> UnitigOcc {
        UnitigOcc::decode_piscem(self.ref_shift, self.pos_mask, word)
    }
}

impl U2Pos for PiscemUnitigTable {
    type EncodedOccs<'a> = OccsRange where Self: 'a;

    fn encoded_unitig_occs(&self, ui: usize) -> Self::EncodedOccs<'_> {
        let start_idx = self.contig_offsets.get(ui) as usize;
        let end_idx = self.contig_offsets.get(ui + 1) as usize;
        // let r = start_idx..end_idx;
        OccsRange::new(start_idx, end_idx)
    }

    fn decode_unitig_occs<T>(&self, _k2u: &T, occs: Self::EncodedOccs<'_>) -> Vec<UnitigOcc> {
        occs.as_range()
            .map(|i| {
                let word = self.ctable.get(i);
                self.decode_occ_word(word)
            })
            .collect()
    }

    fn n_unitig_occs(&self, occs: Self::EncodedOccs<'_>) -> usize {
        occs.len()
    }

    fn empty_encoded_occs(&self) -> Self::EncodedOccs<'_> {
        OccsRange::empty()
    }
}
