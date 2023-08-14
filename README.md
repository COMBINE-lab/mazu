# mazu
A Rust library for building modular, fast and compact indexes over genomic data

> _Mazu ([媽祖](https://en.wikipedia.org/wiki/Mazu))... revered as a tutelary deity of seafarers, including fishermen and sailors..._

## Disclaimer --- This library is in _alpha_ and is under active development.

## Highlights
1. Query ready indexes via plug-and-play k-mer-to-unitig and unitig-to-occurrence mappings.
2. Load (only) compatibility with [pufferfish](https://github.com/COMBINE-lab/pufferfish), deserialize pufferfish indices and work with them in Rust.
3. Streaming queries for generic indexes for free with `.as_streaming()`
4. An easy test-bed for new compression algorithms for unitig-occurrences and k-mer dictionaries.
5. No more CMake.

## Examples

```Rust
// Load a pufferfish index from C++ implementation
let p = to_abs_path(YEAST_CHR01_INDEX);
let pi = DenseIndex::deserialize_from_cpp(p).unwrap();
// Extract unitigs and build a SSHash
let unitig_set = pi.as_ref().clone();
let sshash = SSHash::from_unitig_set(unitig_set, 15, 32, WyHashState::default()).unwrap();

// Drop in an SSHash for a new index
let pi = ModIndex::from_parts(
    pi.base.clone(),
    sshash,
    pi.as_u2pos().clone(),
    pi.as_refseqs().clone(),
);

// Generic implementations take care of query and validation
pi.validate_self();

// Attach a streaming cache and drive the index.
let driver = pi.as_streaming();
driver.validate_self();
```
