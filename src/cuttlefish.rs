use serde::{Deserialize, Serialize};
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    str::FromStr,
};

use crate::Result;

/// Extension for cuttlefish _segments_ file containing unitigs
const CF_SEGS_EXT: &str = "cf_seg";

/// Extension for cuttlefish _sequence_ file containing tiling of reference sequences
const CF_SEQS_EXT: &str = "cf_seq";

// Files for cuttlefish GFA reduced format
#[derive(Debug)]
pub struct CfFiles {
    pub segs: PathBuf,
    pub tiling: PathBuf,
    pub json: PathBuf,
}

// Example .json from cuttlefish
// {
//     "parameters info": {
//         "input": "GCA_000001405.28_GRCh38.p13_genomic.fna.gz",
//         "k": 31,
//         "output prefix": "grch38"
//     },
//     "basic info": {
//         "vertex count": 2505448692
//     },
//     "contigs info": {
//         "maximal unitig count": 36145130,
//         "vertex count in the maximal unitigs": 2505448692,
//         "shortest maximal unitig length": 31,
//         "longest maximal unitig length": 33414,
//         "sum maximal unitig length": 3589802592,
//         "avg. maximal unitig length": 99,
//         "_comment": "lengths are in bases"
//     }
// }
#[derive(Debug, Serialize, Deserialize)]
pub struct CfInfo {
    #[serde(alias = "parameters info")]
    pub parameter_info: CfParamInfo,
    #[serde(alias = "basic info")]
    pub basic_info: CfBasicInfo,
    #[serde(alias = "contigs info")]
    pub contigs_info: CfContigsInfo,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CfParamInfo {
    #[serde(alias = "input")]
    pub input: PathBuf,
    #[serde(alias = "k")]
    pub k: usize,
    #[serde(alias = "output prefix")]
    pub output_prefix: PathBuf,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CfContigsInfo {
    #[serde(alias = "maximal unitig count")]
    pub maximal_unitig_count: usize,
    #[serde(alias = "vertex count in the maximal unitigs")]
    pub vertex_count_in_maximal_unitigs: usize,
    #[serde(alias = "shortest maximal unitig length")]
    pub shortest_maximal_unitig_length: usize,
    #[serde(alias = "longest maximal unitig length")]
    pub longest_maximal_unitig_length: usize,
    #[serde(alias = "sum maximal unitig length")]
    pub sum_maximal_unitig_length: usize,
    #[serde(alias = "avg. maximal unitig length")]
    pub avg_maximal_unitig_length: usize,
    #[serde(alias = "_comment")]
    pub comment: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CfBasicInfo {
    #[serde(alias = "vertex count")]
    pub vertex_count: usize,
}

impl CfInfo {
    pub fn from_path<P: AsRef<Path>>(fp: P) -> crate::err::Result<Self> {
        let fp = fp.as_ref();
        assert!(fp.exists());
        let f = std::fs::File::open(fp).expect("Could not open info file");
        Ok(serde_json::from_reader(f)?)
    }

    pub fn k(&self) -> usize {
        self.parameter_info.k
    }

    pub fn n_unitigs(&self) -> usize {
        self.contigs_info.maximal_unitig_count
    }

    pub fn total_len(&self) -> usize {
        self.contigs_info.sum_maximal_unitig_length
    }
}

impl CfFiles {
    pub fn new<P: AsRef<Path>>(prefix: P) -> Self {
        let prefix = prefix.as_ref().to_path_buf();
        let segs = prefix.with_extension(CF_SEGS_EXT);
        let tiling = prefix.with_extension(CF_SEQS_EXT);
        let json = prefix.with_extension("json");

        Self { segs, tiling, json }
    }
}

use crate::Orientation;

#[derive(Clone, Debug, PartialEq)]
pub enum CfSeqToken {
    PolyN(usize),
    Unitig { id: usize, o: Orientation },
}

impl FromStr for CfSeqToken {
    type Err = crate::Error;
    fn from_str(s: &str) -> Result<Self> {
        if let Some(stripped_s) = s.strip_prefix('N') {
            let len = stripped_s
                .parse()
                .map_err(|_| crate::Error::CfSeqTokenParseError)?;
            Ok(Self::PolyN(len))
        } else {
            let end = s.len() - 1;
            let o = if s.as_bytes()[end] == b'+' {
                Orientation::Forward
            } else {
                Orientation::Backward
            };

            let id = s[..end]
                .parse()
                .map_err(|_| crate::Error::CfSeqTokenParseError)?;
            Ok(Self::Unitig { id, o })
        }
    }
}

// Iterator for cuttlefish reduced GFA format
// Returns Vec<CfSeqTokens> for each line in .seq file.
// Vec of tokens follows reduced GFA format where unitigs or polyN stretches are space separated
// which are:
// - <cfid><o>, where cfid is unique ID cuttlefish assigns to a unitig, and o = + if forward and - otherwise
// - N<n>, for a polyN stretch of length n
pub struct CfSeqIterator<P> {
    lines: std::io::Lines<P>,
}

impl<P: BufRead> Iterator for CfSeqIterator<P> {
    type Item = (String, Vec<CfSeqToken>);

    fn next(&mut self) -> Option<Self::Item> {
        let line = self.lines.next()?;
        let line = line.unwrap();
        let (id, tokens) = line.split_once('\t').unwrap();

        let tokens = tokens.split(' ').map(|s| s.parse().unwrap()).collect();
        Some((id.to_string(), tokens))
    }
}

impl CfFiles {
    pub fn iter_tiling(&self) -> std::io::Result<CfSeqIterator<BufReader<File>>> {
        let f = File::open(self.tiling.clone())?;
        let pb = BufReader::new(f);
        Ok(CfSeqIterator { lines: pb.lines() })
    }
}
