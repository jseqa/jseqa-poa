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
java -jar target/jseqa-poa-0.1.0-SNAPSHOT-jar-with-dependencies.jar \
  --input reads.fasta --type fitting --output alignment.txt --dot graph.dot
```

`--type` accepts `global`, `local`, `fitting`, or `overlap`. Without
`--output`, aligned rows are written one per line to standard output.

To construct the graph directly from a CLUSTAL-format MSA, use `--alignment`
instead of `--input`:

```shell
java -jar target/jseqa-poa-0.1.0-SNAPSHOT-jar-with-dependencies.jar \
  --alignment existing.aln --type fitting --output alignment.txt
```

`--input` and `--alignment` are mutually exclusive.

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
