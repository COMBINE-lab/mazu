use std::fs::File;
use std::path::PathBuf;

use anyhow::bail;
use anyhow::Context;
use anyhow::Result;
use clap::Parser;

const SSHASH_SKEW_DEFAULT: usize = 64;
const PISCEM_EXT: &str = "piscem";
const PUFFERFISH_DENSE_EXT: &str = "pf_dense";

use mazu::index::defaults::{PiscemIndexDefault, PufferfishDenseIndexDefault};
use mazu::Validate;

#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, Clone, Parser)]
enum CLI {
    #[command(subcommand)]
    Build(BuildArgs),

    ValidateFasta(ValidateFastaArgs),
    // Stats(StatsArgs),
}

#[derive(Debug, Clone, Parser)]
enum BuildArgs {
    Piscem(PiscemArgs),
    Pufferfish(PufferfishArgs),
} // Build Piscem or PF1

#[derive(Debug, Clone, Parser)]
struct PiscemArgs {
    #[clap(short, long, default_value_t = 8)]
    n_threads: usize,

    #[clap(short = 'i', long)]
    cuttlefish_prefix: PathBuf,

    #[clap(short, long)]
    output_prefix: PathBuf,

    #[clap(short, long)]
    minimizer_size: usize,
}

#[derive(Debug, Clone, Parser)]
struct PufferfishArgs {
    #[clap(short, long, default_value_t = 8)]
    n_threads: usize,

    #[clap(short = 'i', long)]
    cuttlefish_prefix: PathBuf,

    #[clap(short, long)]
    output_prefix: PathBuf,
}

#[derive(Debug, Clone, Parser)]
struct ValidateFastaArgs {
    #[clap(short, long)]
    index: PathBuf,

    #[clap(short, long)]
    fasta: PathBuf,
}

// #[derive(Debug, Clone, Parser)]
// struct StatsArgs {}

// struct SampleUnitigOccsArgs{} // sample unitig occs for PF1 or PF2 or piscem index.

fn main() -> Result<()> {
    let args = CLI::parse();

    simple_logger::SimpleLogger::new().env().init().unwrap();

    match &args {
        CLI::Build(args) => build(args),
        CLI::ValidateFasta(args) => validate(args),
        // CLI::Stats(args) => todo!(),
    }
}

fn build(args: &BuildArgs) -> Result<()> {
    match args {
        BuildArgs::Piscem(args) => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(args.n_threads)
                .build_global()
                .unwrap();

            log::info!(
                "Building Piscem index from {}.*",
                args.cuttlefish_prefix.display()
            );
            let index = PiscemIndexDefault::from_cf_prefix(
                &args.cuttlefish_prefix,
                args.minimizer_size,
                SSHASH_SKEW_DEFAULT,
            )?;

            let outname = PathBuf::from(&args.output_prefix).with_extension(PISCEM_EXT);
            let outfile = File::create(&outname)?;
            bincode::serialize_into(outfile, &index)?;
            log::info!("Serialized Piscem index to {}", outname.display());
            Ok(())
        }

        BuildArgs::Pufferfish(args) => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(args.n_threads)
                .build_global()
                .unwrap();

            log::info!(
                "Building Pufferfish (Dense) index from {}.*",
                args.cuttlefish_prefix.display()
            );
            let index = PufferfishDenseIndexDefault::from_cf_prefix(&args.cuttlefish_prefix)?;

            let outname = PathBuf::from(&args.output_prefix).with_extension(PUFFERFISH_DENSE_EXT);
            let outfile = File::create(&outname)?;
            bincode::serialize_into(outfile, &index)?;
            log::info!(
                "Serialized Pufferfish (Dense index to {}",
                outname.display()
            );
            Ok(())
        }
    }
}

fn validate(args: &ValidateFastaArgs) -> Result<()> {
    let index_fp = &args.index;
    let fasta = &args.fasta;
    let ext = index_fp
        .extension()
        .with_context(|| "No extension found")?
        .to_str()
        .with_context(|| "Could not parse extension")?;

    match ext {
        PUFFERFISH_DENSE_EXT => {
            let reader = File::open(index_fp)?;
            log::info!(
                "Loading Pufferfish (Dense) index from {}",
                index_fp.display()
            );
            let index: PufferfishDenseIndexDefault = bincode::deserialize_from(reader)?;

            log::info!("Validating against: {}", fasta.display());
            index.validate_fasta(fasta)
        }
        PISCEM_EXT => {
            let reader = File::open(index_fp)?;
            log::info!("Loading Piscem index from {}", index_fp.display());
            let index: PiscemIndexDefault = bincode::deserialize_from(reader)?;
            log::info!("Validating against: {}", fasta.display());
            index.validate_fasta(fasta)
        }
        _ => {
            bail!("Invalid extension: {ext}")
        }
    }
    Ok(())
}
