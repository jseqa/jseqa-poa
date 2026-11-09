# jSeqA-poa

Java partial-order alignment (POA) support for the `jseqa` sequence-alignment
library. The project provides global, local, fitting, and overlap alignment
graphs and can emit a multiple sequence alignment, consensus sequence, or DOT
graph.

## Build and test

The project targets Java 8 or later and uses Maven:

```shell
mvn verify
```

## Command line

Build the executable jar and align sequences from a FASTA file:

```shell
mvn package
java -jar target/jseqa-poa.jar \
  --input reads.fasta --type fitting --output alignment.txt --dot graph.dot
```

`--type` accepts `global`, `local`, `fitting`, or `overlap`. Without
`--output`, aligned rows are written one per line to standard output.

To construct the graph directly from a CLUSTAL-format MSA, use `--alignment`
instead of `--input`:

```shell
java -jar target/jseqa-poa.jar \
  --alignment existing.aln --type fitting --output alignment.txt
```

`--input` and `--alignment` are mutually exclusive.

### Similarity-guided agglomerative POA

The default FASTA workflow remains input-order progressive alignment. To build
a normalized pairwise-similarity guide tree and merge POA profiles bottom-up:

```shell
java -jar target/jseqa-poa.jar \
  --input references.fasta \
  --output references.aln \
  --progressive agglomerative \
  --linkage upgma \
  --pairwise-mode global \
  --score-normalization self-score \
  --threads 4 \
  --guide-tree-output references.nwk
```

The library API returns the merged graph, MSA, guide tree, and the mapping from
MSA rows to original input indices:

```java
AgglomerativePOABuilder builder = new AgglomerativePOABuilder(
        scoreMatrix,
        Linkage.UPGMA,
        ScoreNormalization.SELF_SCORE)
        .pairwiseMode(PairwiseAlignmentMode.GLOBAL)
        .threads(4);

ProgressivePOAResult result = builder.build(sequences);
PartialOrderAlignment graph = result.graph();
List<String> msa = result.alignment();
String newick = result.guideTree().toNewick();
List<Integer> rowToInput = result.alignmentInputIndices();
```

Agglomerative construction aligns complete child profiles with affine gaps; it
does not flatten a subtree into sequences that are re-added one by one.
Duplicate sequences remain distinct observations. Pairwise scoring can use
global, overlap, fitting, or local boundaries, while profile merging is global.

### Seed-derived DAG banding

`SeededBandedPOAAligner` provides an independent read-only global or local
aligner that uses exact query/observed-graph-path seeds to build a compact
corridor through the DAG. Existing dense graph APIs are unchanged:

```java
SeededBandOptions options = SeededBandOptions.builder()
        .kmerLength(11)
        .seedSelection(SeedSelection.MINIMIZERS)
        .minimizerWindow(8)
        .queryPadding(24)
        .graphPadding(8)
        .fallback(SeededBandFallback.DENSE_IF_UNREACHABLE)
        .build();

SeededBandedPOAAligner aligner = new SeededBandedPOAAligner(
        graph, scoreMatrix, SeededGraphMode.IMMUTABLE);
SeededPOAAlignmentResult result = aligner.align(query, options);
List<POGNode> path = result.path();
SeededBandStats diagnostics = result.stats();
```

Use `SeededGraphMode.IMMUTABLE` when aligning many queries to a fixed graph. It
constructs the topology once and skips the otherwise per-query graph fingerprint
scan. The default `MUTABLE` mode detects graph structure or node/read membership
changes and invalidates cached topology and seed indexes automatically. Modifying
a graph after constructing an `IMMUTABLE` aligner is unsupported.

Seed indexes are cached by their k-mer and selection settings. Seed chaining
requires DAG reachability, and the compact
workspace stores only row intervals inside the derived corridor. Since a
seed-derived band is heuristic, fallback policy is explicit in the options.
For `POGraphLocal`, restart cells, best-cell endpoint selection, traceback, and
the returned query interval follow the existing dense local semantics.

A fixture benchmark reports throughput, evaluated/dense cells, and fallbacks:

```shell
mvn package test-compile
java -cp target/test-classes:target/jseqa-poa.jar \
  jseqa.poa.cli.SeededBandedPOABenchmark reference.aln queries.fasta 11 24 local
```

The library rejects null and empty sequences. Alignment matrices are sized to
the current graph and query, so practical limits are determined by available
memory.

## Retaining read positions

By default, graph nodes retain the indices of reads assigned to them. To also
retain the zero-based position within each original read, construct any graph
type with a `PositionedPOGNodeFactory`:

```java
POGraphFitting graph = new POGraphFitting(
        "ACGT",
        new PositionedPOGNodeFactory());
graph.add("ACAGT");

for (POGNode node : graph.getGraph().vertexSet()) {
    PositionedPOGNode positioned = (PositionedPOGNode) node;
    System.out.println(positioned.getCH() + " " + positioned.readPositions());
}
```

Each `PositionedPOGNode` maps a read index to exactly one original read
position. Repeating the same mapping is harmless; attempting to associate the
same read with a different position in that node is rejected.

## Building a graph from an existing MSA

All four graph types accept a collection of equal-width aligned rows. A hyphen
(`-`) represents a gap, and the supplied alignment columns are imported
directly rather than realigned:

```java
List<String> msa = Arrays.asList(
        "AC-GT",
        "ACTGT",
        "AC-GC");

POGraphFitting graph = new POGraphFitting(
        msa,
        new PositionedPOGNodeFactory());
```

Equal residues in one column share a graph node; different residues in that
column become aligned alternatives. Edges follow each read across gap runs.
For positioned nodes, coordinates are zero-based positions in the original
ungapped read. Null rows, unequal row widths, empty/all-gap rows, and all-gap
columns are rejected.

## Performance benchmark

The test sources include a repeatable single-thread benchmark for a human IGHV
reference POA and 350–400 nucleotide queries. It warms up the JVM, excludes
file loading and graph construction from measurement, checks path and query
interval stability, and reports throughput, allocation per read, and GC count:

```shell
mvn package test-compile
java -cp target/test-classes:target/jseqa-poa.jar \
  jseqa.poa.cli.POAPerformanceBenchmark \
  /path/to/V.aln /path/to/queries.fna 5
```
