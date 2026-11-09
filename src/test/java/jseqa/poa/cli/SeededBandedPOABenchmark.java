package jseqa.poa.cli;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import jseqa.alignment.substitution.CustomDnaSubstitutionMatrix;
import jseqa.fastq.Entry;
import jseqa.fastq.FastaReader;
import jseqa.poa.AbstractPartialOrderGraph;
import jseqa.poa.POGraphGlobal;
import jseqa.poa.POGraphLocal;
import jseqa.poa.POGNode;
import jseqa.poa.SeedSelection;
import jseqa.poa.SeededBandFallback;
import jseqa.poa.SeededBandOptions;
import jseqa.poa.SeededBandedPOAAligner;
import jseqa.poa.SeededGraphMode;
import jseqa.poa.SeededPOAAlignmentResult;

/** Repeatable seeded-corridor benchmark using an MSA graph and FASTA queries. */
public final class SeededBandedPOABenchmark {
    public static void main(String[] args) throws Exception {
        if (args.length < 2) {
            System.err.println("usage: SeededBandedPOABenchmark reference.aln queries.fasta "
                    + "[k] [padding] [global|local] [graph-padding]");
            return;
        }
        Path reference = Paths.get(args[0]);
        Path queriesPath = Paths.get(args[1]);
        int k = args.length > 2 ? Integer.parseInt(args[2]) : 11;
        int padding = args.length > 3 ? Integer.parseInt(args[3]) : 24;
        String mode = args.length > 4 ? args[4] : "global";
        int graphPadding = args.length > 5 ? Integer.parseInt(args[5]) : 8;
        List<String> queries = readFasta(queriesPath);
        CustomDnaSubstitutionMatrix scores =
                new CustomDnaSubstitutionMatrix(5, -4, -4, -1);
        AbstractPartialOrderGraph graph = mode.equalsIgnoreCase("local")
                ? new POGraphLocal(POAlign.readClustalAlignment(reference), scores)
                : new POGraphGlobal(POAlign.readClustalAlignment(reference), scores);
        SeededBandedPOAAligner aligner = new SeededBandedPOAAligner(graph, scores,
                SeededGraphMode.IMMUTABLE);
        SeededBandOptions options = SeededBandOptions.builder()
                .kmerLength(k).seedSelection(SeedSelection.MINIMIZERS)
                .minimizerWindow(8).queryPadding(padding).graphPadding(graphPadding)
                .fallback(SeededBandFallback.DENSE_IF_UNREACHABLE).build();

        long evaluated = 0;
        long dense = 0;
        int fallbacks = 0;
        int noChains = 0;
        List<String> seededPaths = new ArrayList<String>(queries.size());
        long started = System.nanoTime();
        for (String query : queries) {
            SeededPOAAlignmentResult result = aligner.align(query, options);
            seededPaths.add(pathSignature(result.path()));
            evaluated += result.stats().evaluatedCells();
            dense += result.stats().denseCellCount();
            if (result.stats().denseFallbackUsed()) fallbacks++;
            if (result.stats().noUsableChain()) noChains++;
        }
        double seconds = (System.nanoTime() - started) / 1_000_000_000.0;
        int denseAgreements = 0;
        long denseStarted = System.nanoTime();
        for (int index = 0; index < queries.size(); index++) {
            String densePath = pathSignature(graph.alignToGraph(queries.get(index)));
            if (seededPaths.get(index).equals(densePath)) {
                denseAgreements++;
            }
        }
        double denseSeconds = (System.nanoTime() - denseStarted) / 1_000_000_000.0;
        System.out.printf("mode: %s, graph nodes: %d, queries: %d, k: %d, padding: %d, "
                        + "graph padding: %d%n",
                mode, graph.numNodes(), queries.size(), k, padding, graphPadding);
        System.out.printf("wall time: %.3f s, throughput: %.1f reads/s%n",
                seconds, queries.size() / seconds);
        System.out.printf("evaluated cells: %d / %d (%.2f%%)%n", evaluated, dense,
                dense == 0 ? 0.0 : evaluated * 100.0 / dense);
        System.out.printf("dense fallbacks: %d, no-chain queries: %d%n",
                fallbacks, noChains);
        System.out.printf("dense wall time: %.3f s, band/dense path agreements: %d/%d%n",
                denseSeconds, denseAgreements, queries.size());
    }

    private static List<String> readFasta(Path path) throws Exception {
        List<String> sequences = new ArrayList<String>();
        try (FastaReader reader = new FastaReader(path)) {
            for (Entry entry : reader) sequences.add(entry.sequence());
        }
        return sequences;
    }

    private static String pathSignature(List<POGNode> path) {
        StringBuilder signature = new StringBuilder(path.size());
        for (POGNode node : path) signature.append(node == null ? '-' : node.getCH());
        return signature.toString();
    }

    private SeededBandedPOABenchmark() {}
}
