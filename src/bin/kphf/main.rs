use clap::Parser;
use kmers::naive_impl::CanonicalKmerIterator;
use mazu::{
    cuttlefish::*,
    kphf::{pfhash::PFHashDefault, sshash::SSHashDefault, PFHash, SSHash, WyHashState, K2U},
    unitig_set::UnitigSet,
    util::FastaReader,
};
use std::{
    fs::File,
    io::{BufReader, BufWriter},
    path::PathBuf,
};

#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, Clone, Parser)]
pub struct CLI {
    #[clap(short, long, default_value_t = 8)]
    n_threads: usize,

    #[clap(subcommand)]
    cmd: Cmd,
}

#[derive(Debug, Clone, Parser)]
enum Cmd {
    Build(BuildArgs),
    Validate(ValidateArgs),
    Stats(StatsArgs),
    Bench(BenchArgs),
}

#[derive(Debug, Clone, Parser)]
struct StatsArgs {
    #[clap(short, long)]
    input: PathBuf,
}

#[derive(Debug, Clone, Parser)]
struct BenchArgs {
    #[clap(short, long)]
    index: PathBuf,

    #[clap(short, long)]
    fasta: PathBuf,

    #[clap(short, long, default_value_t = false)]
    use_streaming: bool,
}

#[derive(Debug, Clone, Parser)]
struct ValidateArgs {
    #[clap(short, long)]
    input: PathBuf,
    // #[clap(short, long)]
    // use_streaming: bool, // TODO streaming validation
}

#[derive(Debug, Clone, Parser)]
struct BuildArgs {
    #[clap(subcommand)]
    hash: BuildHash,
}

#[derive(Debug, Clone, Parser)]
enum BuildHash {
    SSHash {
        #[clap(short, long)]
        input_prefix: PathBuf,

        #[clap(short, long)]
        output_prefix: PathBuf,

        #[clap(short, long)]
        minimizer_size: usize,

        #[clap(short, long)]
        skew_min: usize,

        #[clap(short, long)]
        validate: bool,
    },

    PFHash {
        #[clap(short, long)]
        input_prefix: PathBuf,

        #[clap(short, long)]
        output_prefix: PathBuf,

        #[clap(short, long)]
        validate: bool,
    },
}

const SSHASH_EXT: &str = "sshash";
const PFHASH_EXT: &str = "pfhash";

use anyhow::{bail, Context, Result};

fn main() -> Result<()> {
    let args = CLI::parse();
    simple_logger::SimpleLogger::new().env().init().unwrap();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.n_threads)
        .build_global()
        .unwrap();

    log::info!(
        "{} Rayon threads in thread pool",
        rayon::current_num_threads()
    );

    match args.cmd {
        Cmd::Build(args) => build(args),
        Cmd::Validate(args) => validate(args),
        Cmd::Stats(args) => stats(args),
        Cmd::Bench(args) => bench(args),
    }
}

fn build(args: BuildArgs) -> Result<()> {
    match args.hash {
        BuildHash::PFHash {
            input_prefix,
            output_prefix,
            validate,
        } => {
            let cf_files = CfFiles::new(input_prefix);

            log::info!("*** Building Unitig Set");
            let (unitigs, _uid_to_idx) = UnitigSet::from_cf_reduced_gfa(&cf_files).unwrap();

            let out_fp = output_prefix.with_extension(PFHASH_EXT);
            let outfile = File::create(&out_fp)?;
            let outfile = BufWriter::new(outfile);

            log::info!("*** Building Pufferfish Hash");

            let hash = PFHash::from_unitig_set(unitigs);
            if validate {
                log::info!("*** Validating");
                let start = std::time::Instant::now();
                hash.validate_self_parallel();
                let elapsed = start.elapsed();
                log::info!("\t {}s to validate", elapsed.as_secs_f64());
                let nanos = elapsed.as_nanos() as f64;
                let ns_per_km = nanos / hash.n_kmers() as f64;
                log::info!("\t {}ns per km to validate", ns_per_km);

                log::info!("\t OK!");
            }

            log::info!("*** Serializing to {}", out_fp.to_string_lossy());
            bincode::serialize_into(outfile, &hash)?;
        }

        BuildHash::SSHash {
            input_prefix,
            output_prefix,
            minimizer_size,
            skew_min,
            validate,
        } => {
            let cf_files = CfFiles::new(input_prefix);

            log::info!("*** Building Unitig Set");
            let (unitigs, _uid_to_idx) = UnitigSet::from_cf_reduced_gfa(&cf_files).unwrap();

            let out_fp = output_prefix.with_extension(SSHASH_EXT);
            let outfile = File::create(&out_fp)?;
            let outfile = BufWriter::new(outfile);

            log::info!("*** Building SSHash");

            let bh = WyHashState::with_seed(0);
            let hash = SSHash::from_unitig_set(unitigs, minimizer_size, skew_min, bh)?;
            if validate {
                log::info!("*** Validating");
                let start = std::time::Instant::now();
                hash.validate_self_parallel();
                let elapsed = start.elapsed();
                log::info!("\t {}s to validate", elapsed.as_secs_f64());
                let nanos = elapsed.as_nanos() as f64;
                let ns_per_km = nanos / hash.n_kmers() as f64;
                log::info!("\t {}ns per km to validate", ns_per_km);
                log::info!("\t OK!");
            }
            hash.print_stats();

            log::info!("*** Serializing to {}", out_fp.to_string_lossy());
            bincode::serialize_into(outfile, &hash)?;
        }
    };
    log::info!("DONE.");

    Ok(())
}

fn validate(args: ValidateArgs) -> Result<()> {
    let input = args.input;
    let ext = input
        .extension()
        .with_context(|| format!("Could not retrieve file extension from {}", input.display()))?;

    match ext.to_str().context("Canot match ext")? {
        SSHASH_EXT => {
            let fp = File::open(input)?;
            log::info!("Loading SSHash");
            let hash: SSHashDefault = bincode::deserialize_from(fp)?;
            log::info!("Loaded!");

            log::info!("Starting validation...");
            let start = std::time::Instant::now();
            hash.validate_self_parallel();
            let elapsed = start.elapsed();
            log::info!("\t {}s to validate", elapsed.as_secs_f64());
            let nanos = elapsed.as_nanos() as f64;
            let ns_per_km = nanos / hash.n_kmers() as f64;
            log::info!("\t {}ns per km to validate", ns_per_km);
            log::info!("\t OK!");
            Ok(())
        }
        PFHASH_EXT => {
            let fp = File::open(input)?;
            log::info!("Loading PFHash");
            let hash: PFHashDefault = bincode::deserialize_from(fp)?;
            log::info!("Loaded!");

            log::info!("Starting validation...");
            let start = std::time::Instant::now();
            hash.validate_self_parallel();
            let elapsed = start.elapsed();
            log::info!("\t {}s to validate", elapsed.as_secs_f64());
            let nanos = elapsed.as_nanos() as f64;
            let ns_per_km = nanos / hash.n_kmers() as f64;
            log::info!("\t {}ns per km to validate", ns_per_km);
            log::info!("\t OK!");
            Ok(())
        }

        _ => bail!("{:?} does not have compatible extension.", input),
    }
}

fn stats(args: StatsArgs) -> Result<()> {
    let input = args.input;
    let ext = input
        .extension()
        .with_context(|| format!("Could not retrieve file extension from {}", input.display()))?;

    match ext.to_str().context("Canot match ext")? {
        SSHASH_EXT => {
            let fp = File::open(input)?;
            let sshash: SSHashDefault = bincode::deserialize_from(fp)?;
            sshash.print_stats();
            Ok(())
        }
        PFHASH_EXT => {
            // let fp = File::open(input)?;
            // let hash: PFHashDefault = bincode::deserialize_from(fp)?;
            // hash.print_stats();
            // Ok(())

            todo!("Not supported yet")
        }

        _ => bail!("{:?} does not have compatible extension.", input),
    }
}

fn bench(args: BenchArgs) -> Result<()> {
    let index = args.index;
    let fasta = args.fasta;
    let ext = index
        .extension()
        .with_context(|| format!("Could not retrieve file extension from {}", index.display()))?;

    log::info!("Starting single threaded benchmark");

    let mut n_kmers = 0_usize;
    let mut n_hit = 0_usize;
    let mut n_miss = 0_usize;
    let start;

    match ext.to_str().context("Canot match ext")? {
        SSHASH_EXT => {
            let fp = File::open(index)?;
            log::info!("Loading SSHash");
            let sshash: SSHashDefault = bincode::deserialize_from(fp)?;
            log::info!("Loaded!");

            let f = File::open(fasta)?;
            let f = BufReader::new(f);
            let rdr = FastaReader::new(f);
            start = std::time::Instant::now();

            if args.use_streaming {
                let mut sshash = sshash.as_streaming();
                rdr.into_iter().for_each(|rec| {
                    let seq = rec.seq.as_bytes();
                    for tup in CanonicalKmerIterator::from_u8_slice(seq, sshash.k() as u8) {
                        n_kmers += 1;
                        match sshash.k2u_streaming(&tup.km) {
                            Some(_) => n_hit += 1,
                            None => n_miss += 1,
                        }
                    }
                });
            } else {
                rdr.into_iter().for_each(|rec| {
                    let seq = rec.seq.as_bytes();
                    for tup in CanonicalKmerIterator::from_u8_slice(seq, sshash.k() as u8) {
                        n_kmers += 1;
                        match sshash.k2u(&tup.km) {
                            Some(_) => n_hit += 1,
                            None => n_miss += 1,
                        }
                    }
                });
            }
        }
        PFHASH_EXT => {
            todo!("Not supported yet");
        }

        _ => bail!("{:?} does not have compatible extension.", index),
    }

    let elapsed = start.elapsed();
    log::info!("\t Completed in {}s", elapsed.as_secs_f64());
    let nanos = elapsed.as_nanos() as f64;
    let ns_per_km = nanos / n_kmers as f64;
    log::info!("\t found {} / {} kmers", n_hit, n_kmers);
    log::info!("\t {}ns per km", ns_per_km);

    Ok(())
}
