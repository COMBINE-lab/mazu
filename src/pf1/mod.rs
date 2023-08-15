mod boophf;
pub mod cpp;
pub mod dense_index;
pub mod sparse_index;
pub(crate) mod test_utils;
mod unitig_table;

use std::{
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
};

use log::debug;

use kmers::naive_impl::seq_vector::SeqVector;
use serde::{Deserialize, Serialize};

use super::{index::BaseIndex, refseq::RefSeqCollection};

use crate::{IndexMetadata, ModIndexType, Result};
use cpp::{DeserializeFromCpp, FromCereal, FromCompact};

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone, Copy)]
pub enum PF1Type {
    #[serde(alias = "dense")] // aliased to lowercase 'd' for c++ compatibility
    Dense,
    #[serde(alias = "sparse")]
    Sparse,
}

impl From<PF1Type> for ModIndexType {
    fn from(pf1_t: PF1Type) -> Self {
        match pf1_t {
            PF1Type::Dense => Self::PF1Dense,
            PF1Type::Sparse => Self::PF1Sparse,
        }
    }
}

#[derive(Debug)]
pub struct PF1FilePaths {
    // TODO: Absolute paths to all files?

    // Parent directory
    pub prefix: PathBuf,

    // binary format files
    pub complete_ref_lens: PathBuf,
    pub ctable: PathBuf,
    pub ctg_offsets: PathBuf,
    pub mphf: PathBuf, // Serialized BooPHF
    pub pos: PathBuf,
    pub sample_pos: PathBuf,        // [sparse index] sampled positions
    pub presence: PathBuf,          // [sparse index]
    pub extension_lengths: PathBuf, // [sparse inndex]
    pub extension_bases: PathBuf,   // [sparse index]
    pub direction: PathBuf,         // [sparse index]
    pub canonical: PathBuf,         // [sparse index]
    pub rank: PathBuf,              // contig boundaries
    pub ref_accum_lens: PathBuf,
    pub ref_lens: PathBuf,
    pub ref_seq: PathBuf,
    pub seq: PathBuf,

    // Other files
    pub duplicate_clusters_tsv: PathBuf,
    pub info_json: PathBuf,
    pub ref_indexing_log: PathBuf,
}

impl PF1FilePaths {
    pub fn new<P: AsRef<Path>>(dir: P) -> Self {
        let dir = dir.as_ref();
        Self {
            prefix: PathBuf::from(dir),
            complete_ref_lens: dir.join(consts::fp::COMPLETE_REF_LENS),
            ctable: dir.join(consts::fp::CTABLE),
            ctg_offsets: dir.join(consts::fp::CTG_OFFSETS),
            mphf: dir.join(consts::fp::MPHF),
            pos: dir.join(consts::fp::POS),
            sample_pos: dir.join(consts::fp::SAMPLE_POS),
            presence: dir.join(consts::fp::PRESENCE),
            extension_lengths: dir.join(consts::fp::EXTENSION_LENGTHS),
            extension_bases: dir.join(consts::fp::EXTENSION_BASES),
            direction: dir.join(consts::fp::DIRECTION),
            canonical: dir.join(consts::fp::CANONICAL),
            rank: dir.join(consts::fp::RANK),
            ref_accum_lens: dir.join(consts::fp::REF_ACCUM_LENS),
            ref_lens: dir.join(consts::fp::REF_LENS),
            ref_seq: dir.join(consts::fp::REF_SEQ),
            seq: dir.join(consts::fp::SEQ),
            duplicate_clusters_tsv: dir.join(consts::fp::DUPLICATE_CLUSTERS_TSV),
            info_json: dir.join(consts::fp::INFO_JSON),
            ref_indexing_log: dir.join(consts::fp::REF_INDEXING_LOG),
        }
    }

    // pub fn from_dir<P: AsRef<Path>>(dname: P) -> Self {
    //     let dir = dname.as_ref();
    //     // directory check is not checkable
    //     if !dir.exists() {
    //         panic!("The directory {} does not exist", dir.to_str().unwrap());
    //     }

    //     if !dir.is_dir() {
    //         panic!(
    //             "The path {} is not a valid directory",
    //             dir.to_str().unwrap()
    //         )
    //     }

    //     Self::with_dir(dname)
    // }
}
pub mod consts {
    pub mod fp {
        pub const COMPLETE_REF_LENS: &str = "complete_ref_lens.bin";
        pub const CTABLE: &str = "ctable.bin";
        pub const CTG_OFFSETS: &str = "ctg_offsets.bin";
        pub const DUPLICATE_CLUSTERS_TSV: &str = "duplicate_clusters.tsv";
        pub const INFO_JSON: &str = "info.json";
        pub const MPHF: &str = "mphf.bin";
        pub const POS: &str = "pos.bin";
        pub const SAMPLE_POS: &str = "sample_pos.bin";
        pub const PRESENCE: &str = "presence.bin";
        pub const EXTENSION_LENGTHS: &str = "extensionSize.bin";
        pub const EXTENSION_BASES: &str = "extension.bin";
        pub const DIRECTION: &str = "direction.bin";
        pub const CANONICAL: &str = "canonical.bin";
        pub const RANK: &str = "rank.bin";
        pub const REF_ACCUM_LENS: &str = "refAccumLengths.bin";
        pub const REF_INDEXING_LOG: &str = "ref_indexing.log";
        pub const REF_LENS: &str = "reflengths.bin";
        pub const REF_SEQ: &str = "refseq.bin";
        pub const SEQ: &str = "seq.bin";
    }
}

// FIXME: Split info to a legacy LegacyInfo or V1Info that contains info-fields
// that are have no Options. Fields like sample_size and extension_size should
// wholly be encoded in 'sampling_type'
// and 'sampling_type' should be renamed to pufferfish_type
#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct PF1Info {
    pub index_version: u64,
    pub reference_gfa: Vec<String>,
    pub sampling_type: PF1Type,

    #[serde(rename = "k")]
    pub kmer_size: usize,
    pub num_kmers: usize,
    pub num_contigs: usize,
    #[serde(rename = "seq_length")]
    pub seq_len: usize,
    pub have_ref_seq: bool,
    pub have_edge_vec: bool,

    #[serde(rename = "SeqHash")]
    pub seq_hash: String,
    #[serde(rename = "NameHash")]
    pub name_hash: String,
    #[serde(rename = "SeqHash512")]
    pub seq_hash_512: String,
    #[serde(rename = "NameHash512")]
    pub name_hash_512: String,

    #[serde(rename = "DecoySeqHash")]
    pub decoy_seq_hash: String,
    #[serde(rename = "DecoyNameHash")]
    pub decoy_name_hash: String,
    pub num_decoys: usize,

    // only present if this is a sparse index
    pub sample_size: Option<usize>,
    pub extension_size: Option<usize>,

    pub first_decoy_index: usize,
    pub keep_duplicates: bool,
}

impl PF1Info {
    pub fn load<P: AsRef<Path>>(p: P) -> Result<Self> {
        let p = p.as_ref();
        dbg!(p);
        debug!("Loading Info from: {:?}", &p);
        let f = File::open(p)?;

        let r = BufReader::new(f);

        Ok(serde_json::from_reader(r)?)
    }

    pub fn as_base_index(&self) -> BaseIndex {
        let metadata = IndexMetadata {
            have_edge_vec: self.have_edge_vec,
            seq_hash: self.seq_hash.clone(),
            name_hash: self.name_hash.clone(),
            seq_hash_512: self.seq_hash_512.clone(),
            name_hash_512: self.name_hash_512.clone(),
            decoy_seq_hash: self.decoy_seq_hash.clone(),
            decoy_name_hash: self.decoy_name_hash.clone(),
            num_decoys: self.num_decoys,
            first_decoy_index: self.first_decoy_index,
            keep_duplicates: self.keep_duplicates,
        };

        BaseIndex::new()
            .set_index_type(ModIndexType::from(self.sampling_type))
            .set_metadata(metadata)
    }
}

impl DeserializeFromCpp for RefSeqCollection {
    fn deserialize_from_cpp<P: AsRef<Path>>(dir: P) -> Result<Self> {
        let files = PF1FilePaths::new(&dir);

        debug!("Loading refseq");
        let seq = if files.ref_seq.exists() {
            Some(SeqVector::from_compact_serialized(files.ref_seq)?)
        } else {
            None
        };

        debug!("Loading reflens");
        // let ref_lens = Vec::load_from_cereal_archive(files.ref_lens)?;
        let prefix_sum: Vec<u64> = {
            // Note: appending zero to make a valid prefix sum of lengths
            let mut v = Vec::load_from_cereal_archive(files.ref_accum_lens)?;
            v.insert(0, 0);
            v
        };
        // let complete_ref_lens = Vec::load_from_cereal_archive(files.complete_ref_lens)?;
        let prefix_sum = prefix_sum.into_iter().map(|x| x as usize).collect();
        Ok(Self::from_parts(seq, prefix_sum))
    }
}
