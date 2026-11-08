# Changelog

Notable changes to jSeqA-poa are documented here.

## [0.1.0]

### Added

- Independent seed-derived global and local sequence-to-POA alignment with observed-path
  k-mer/minimizer indexing, reachability-aware anchor chaining, compact DAG
  corridor storage, explicit dense fallbacks, and detailed diagnostics.
- Similarity-guided agglomerative POA construction with normalized affine-gap
  pairwise scoring, deterministic UPGMA/single/complete/average guide trees,
  and true profile-to-profile POA merging.
  mappings, including distinct mappings for duplicate sequences.
- Bounded parallel pairwise scoring with deterministic results across thread
  counts.
- Direct construction of global, local, fitting, and overlap POA graphs from
  an existing multiple sequence alignment.
- Optional `PositionedPOGNode` nodes that retain one zero-based position per
  read at each graph node.
- Extensible node creation through `POGNodeFactory`.
- CLI support for CLUSTAL `.aln` input through `--alignment`.

### Fixed

- Retained locally unaligned read prefixes and suffixes in MSA output without
  changing the POA graph, alignment path, or consensus semantics.
- Rejected cyclic alignment-column graphs and duplicate read membership within
  a single output column instead of silently emitting a corrupted alignment.
- Included the final graph row and query column when selecting a local
  alignment endpoint.
- Corrected global start/end-node indexing and handling of single-node reads.
- Corrected source-node mismatch classification, origin pointers, and leading
  query-gap scoring.
- Replaced fixed-size alignment matrices with matrices sized for the current
  graph and query.
- Made graph copies independent, including nodes, aligned-node relationships,
  edge read metadata, and positioned-node metadata.
- Corrected mismatch-node relationships and edge/read bookkeeping during
  backtracking.
- Fixed alignment-column iteration when only aligned alternatives remained.
- Counted source-node support when calculating the consensus path.

### Changed

- Updated Log4j to 2.26.1 and JUnit to 4.13.2.
- Improved CLI validation, error messages, plain-text MSA output, and
  documentation.
- Reduced alignment-loop allocation by using primitive score and pointer
  storage and allocation-free candidate selection.
- Reused grow-only dynamic-programming matrices across queries while retaining
  the active matrix dimensions for each alignment.
- Cached graph topology and primitive predecessor rows until the graph changes.
- Added an immutable seeded-aligner mode that constructs topology once and
  eliminates per-read graph fingerprint traversal and its allocations.
- Added performance-compatibility tests and a repeatable IgPOA fixture
  benchmark that reports throughput, allocation, and garbage collection.

## [0.0.1]

- Initial Java POA implementation built on jSeqA and JGraphT.
- Supported progressive global, local, fitting, and overlap alignment of
  unaligned sequences.
- Provided consensus generation, MSA generation, DOT export, and a FASTA-based
  command-line entry point.
- Stored read membership on nodes and edges, but did not retain original read
  positions.
- Did not support direct graph construction from an existing MSA or CLUSTAL
  input in the CLI.
- Used fixed `2000 x 500` dynamic-programming matrices and older dependency
  versions.
