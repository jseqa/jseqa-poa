# Changelog

Notable changes to jSeqA-poa are documented here.

## [0.1.0] - Unreleased

### Added

- Direct construction of global, local, fitting, and overlap POA graphs from
  an existing multiple sequence alignment.
- Optional `PositionedPOGNode` nodes that retain one zero-based position per
  read at each graph node.
- Extensible node creation through `POGNodeFactory`.
- CLI support for CLUSTAL `.aln` input through `--alignment`; this option is
  mutually exclusive with unaligned FASTA input through `--input`.
- Regression coverage for alignment boundaries, graph copying, consensus,
  matrix sizing, MSA import, positioned nodes, CLI input/output, and DOT output.

### Fixed

- Included the final graph row and query column when selecting a local
  alignment endpoint.
- Corrected global start/end-node indexing and handling of single-node reads.
- Corrected source-node mismatch classification, origin pointers, and leading
  query-gap scoring.
- Replaced fixed-size alignment matrices with matrices sized for the current
  graph and query.
- Removed ten artificial trailing gap columns from generated MSAs.
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
