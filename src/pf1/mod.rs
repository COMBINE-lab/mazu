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

use super::{
    index::{BaseIndex, PufferfishType},
    refseq::RefSeqCollection,
};

use crate::Result;
use cpp::{DeserializeFromCpp, FromCereal, FromCompact};

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

impl DeserializeFromCpp for BaseIndex {
    fn deserialize_from_cpp<P: AsRef<Path>>(dir: P) -> Result<Self> {
        let files = PF1FilePaths::new(dir);

        let info = PF1Info::load(files.info_json);
        //   let ref_seq = Self::load_ref_seq(&info);
        //   let edge_vec = Self::load_edge_vec(&info);

        // debug!("Loading mphf");
        // let mphf = BooPHF::<u64>::deserialize_from_cpp(files.mphf)?;

        // debug!("Loading seq");
        // let seq = SeqVector::from_compact_serialized(files.seq)?;
        // assert_eq!(seq.len(), info.seq_len);

        // let last_seq_pos = seq.len() - (info.kmer_size as usize);

        // debug!("Loading refseq");
        // let ref_seq = if info.have_ref_seq {
        //     Some(SeqVector::from_compact_serialized(files.ref_seq)?)
        // } else {
        //     None
        // };

        // debug!("Loading reflens");
        // //TODO assert that the len of ref_seq is the sum of lens of reference sequences
        // let ref_lens = Vec::load_from_cereal_archive(files.ref_lens)?;
        // let ref_accum_lens = {
        //     // TODO: note, we are appending one zero to the front so accessing the offset is easy
        //     let mut v = Vec::load_from_cereal_archive(files.ref_accum_lens)?;
        //     v.insert(0, 0);
        //     v
        // };
        // let complete_ref_lens = Vec::load_from_cereal_archive(files.complete_ref_lens)?;

        // let (ref_names, ref_exts, ctable) = Self::deserialize_contig_table(files.ctable)?;

        // debug!("Loading bv");
        // let bv = {
        //     let mut bv = BitVector::from_compact_serialized(files.rank)?;
        //     bv.enable_rank();
        //     bv.enable_select();
        //     bv
        // };
        // assert_eq!(bv.len(), info.seq_len);

        debug!("Loaded base index");

        Ok(Self {
            // seq,
            // last_seq_pos,
            // mphf,
            // bv,
            // ref_seq,
            // ref_lens,
            // ref_accum_lens, // cumalative length of references
            // _complete_ref_lens: complete_ref_lens,

            // Info fields
            index_version: info.index_version,
            reference_gfa: info.reference_gfa,
            k: info.kmer_size,
            num_kmers: info.num_kmers,
            num_contigs: info.num_contigs,
            seq_len: info.seq_len,
            // [omitted] have_ref_seq: info.have_ref_seq,
            // [omitted] num_kmers: usize,
            // [omitted] num_contigs: info.num_contigs,
            // [omitted] seq_len: usize,
            // [omitted] have_ref_seq: info.have_ref_seq,
            have_edge_vec: info.have_edge_vec,
            seq_hash: info.seq_hash,
            name_hash: info.name_hash,
            seq_hash_512: info.seq_hash_512,
            name_hash_512: info.name_hash_512,
            decoy_seq_hash: info.decoy_seq_hash,
            decoy_name_hash: info.decoy_name_hash,
            num_decoys: info.num_decoys,
            first_decoy_index: info.first_decoy_index,
            keep_duplicates: info.keep_duplicates,
        })
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
    pub sampling_type: PufferfishType,

    #[serde(rename = "k")]
    pub kmer_size: usize,
    pub num_kmers: usize,
    pub num_contigs: usize,
    #[serde(rename = "seq_length")]
    pub seq_len: usize,
    pub have_ref_seq: bool,  //FIXME: has_ref_seq?
    pub have_edge_vec: bool, //FIXME: has_edge_vec?

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
    pub fn load<P: AsRef<Path>>(p: P) -> Self {
        let p = p.as_ref();
        debug!("Loading Info from: {:?}", &p);
        let f = File::open(p).unwrap();

        let r = BufReader::new(f);

        serde_json::from_reader(r).unwrap()
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
