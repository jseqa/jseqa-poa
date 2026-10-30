package jseqa.poa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import jseqa.alignment.substitution.CustomDnaSubstitutionMatrix;
import org.junit.Test;

public class SeededBandedPOAAlignerTest {
    private final CustomDnaSubstitutionMatrix scores =
            new CustomDnaSubstitutionMatrix(5, -4, -4, -1);

    @Test
    public void exactSeedsProduceCompactExactMatchCorridor() {
        POGraphGlobal graph = new POGraphGlobal("ACGTTGCAACGT", scores);
        SeededBandedPOAAligner aligner = new SeededBandedPOAAligner(graph, scores);
        SeededPOAAlignmentResult result = aligner.align("ACGTTGCAACGT",
                options(4, 1, SeededBandFallback.NEVER));

        assertEquals("ACGTTGCAACGT", pathBases(result.path()));
        assertEquals(60, result.score());
        assertTrue(result.anchors().size() >= 2);
        assertTrue(result.stats().evaluatedCells() < result.stats().denseCellCount());
        assertFalse(result.stats().denseFallbackUsed());
    }

    @Test
    public void seedCorridorTraversesTheCompatibleGraphBranch() {
        POGraphGlobal graph = new POGraphGlobal(Arrays.asList(
                "ACG-TACG",
                "ATGGTACG",
                "ACGGTTCG"), scores);
        SeededBandedPOAAligner aligner = new SeededBandedPOAAligner(graph, scores);
        SeededPOAAlignmentResult result = aligner.align("ATGGTACG",
                options(3, 2, SeededBandFallback.DENSE_IF_UNREACHABLE));

        assertEquals("ATGGTACG", pathBases(result.path()));
        assertTrue(result.stats().chainedAnchors() >= 2);
    }

    @Test
    public void noSharedSeedsFallsBackToDenseAlignment() {
        POGraphGlobal graph = new POGraphGlobal("AAAAAAAAAAAA", scores);
        SeededBandedPOAAligner aligner = new SeededBandedPOAAligner(graph, scores);
        List<POGNode> dense = graph.alignToGraph("CCCCCCCCCCCC");
        SeededPOAAlignmentResult result = aligner.align("CCCCCCCCCCCC",
                options(5, 2, SeededBandFallback.DENSE_IF_NO_CHAIN));

        assertTrue(result.stats().noUsableChain());
        assertTrue(result.stats().denseFallbackUsed());
        assertEquals(-48, result.score());
        assertEquals(dense, result.path());
    }

    @Test(expected = SeededBandAlignmentException.class)
    public void noFallbackRejectsAQueryWithoutSeeds() {
        POGraphGlobal graph = new POGraphGlobal("AAAAAAAA", scores);
        new SeededBandedPOAAligner(graph, scores).align("CCCCCCCC",
                options(4, 2, SeededBandFallback.NEVER));
    }

    @Test
    public void ambiguousQueryKmersAreSkipped() {
        POGraphGlobal graph = new POGraphGlobal("ACGTACGTACGT", scores);
        SeededPOAAlignmentResult result = new SeededBandedPOAAligner(graph, scores)
                .align("NNNNNNNNNNNN", options(4, 2,
                        SeededBandFallback.DENSE_IF_NO_CHAIN));
        assertEquals(0, result.stats().querySeeds());
        assertTrue(result.stats().denseFallbackUsed());
    }

    @Test
    public void seedIndexIsReusedAndInvalidatedAfterGraphMutation() {
        POGraphGlobal graph = new POGraphGlobal("ACGTTGCA", scores);
        SeededBandedPOAAligner aligner = new SeededBandedPOAAligner(graph, scores);
        SeededBandOptions options = options(3, 2, SeededBandFallback.DENSE_IF_NO_CHAIN);
        aligner.align("ACGTTGCA", options);
        assertEquals(1, aligner.cachedSeedIndexCount());
        aligner.align("ACGTTGCA", options);
        assertEquals(1, aligner.cachedSeedIndexCount());

        graph.add("TTGGAACC");
        SeededPOAAlignmentResult result = aligner.align("TTGGAACC", options);
        assertEquals(1, aligner.cachedSeedIndexCount());
        assertTrue(result.stats().graphHits() > 0);
    }

    @Test
    public void immutableModeDoesNotReadTheGraphForEveryQuery() {
        CountingPOGraphGlobal graph = new CountingPOGraphGlobal("ACGTTGCA", scores);
        SeededBandedPOAAligner aligner = new SeededBandedPOAAligner(graph, scores,
                SeededGraphMode.IMMUTABLE);
        assertEquals(1, graph.graphReads);

        aligner.align("ACGTTGCA", options(3, 2, SeededBandFallback.NEVER));
        aligner.align("ACGTTGCA", options(3, 2, SeededBandFallback.NEVER));
        aligner.align("ACGTTGCA", options(3, 2, SeededBandFallback.NEVER));

        assertEquals(1, graph.graphReads);
    }

    @Test
    public void immutableModeBuildsEachSeedIndexOnlyOnce() {
        CountingPOGraphGlobal graph = new CountingPOGraphGlobal("ACGTTGCA", scores);
        SeededBandedPOAAligner aligner = new SeededBandedPOAAligner(graph, scores,
                SeededGraphMode.IMMUTABLE);
        SeededBandOptions allKmers = options(3, 2, SeededBandFallback.NEVER);

        aligner.align("ACGTTGCA", allKmers);
        aligner.align("ACGTTGCA", allKmers);
        assertEquals(1, aligner.cachedSeedIndexCount());

        SeededBandOptions minimizers = SeededBandOptions.builder()
                .kmerLength(3).seedSelection(SeedSelection.MINIMIZERS)
                .minimizerWindow(2).queryPadding(2).graphPadding(1)
                .minimumChainedSeeds(2).fallback(SeededBandFallback.NEVER).build();
        aligner.align("ACGTTGCA", minimizers);
        aligner.align("ACGTTGCA", minimizers);
        assertEquals(2, aligner.cachedSeedIndexCount());
        assertEquals(1, graph.graphReads);
    }

    @Test(expected = SeededBandAlignmentException.class)
    public void immutableModeDeliberatelyDoesNotObserveLaterGraphMutations() {
        POGraphGlobal graph = new POGraphGlobal("AAAAAAAAAAAA", scores);
        SeededBandedPOAAligner aligner = new SeededBandedPOAAligner(graph, scores,
                SeededGraphMode.IMMUTABLE);
        graph.add("CCCCCCCCCCCC");

        aligner.align("CCCCCCCCCCCC",
                options(4, 2, SeededBandFallback.NEVER));
    }

    @Test
    public void repeatedSeedsCanBeFilteredAndTriggerFallback() {
        POGraphGlobal graph = new POGraphGlobal(Arrays.asList(
                "AAAAAAAAAA", "AAAAAAAAAA", "AAAAAAAAAA"), scores);
        SeededBandOptions options = SeededBandOptions.builder()
                .kmerLength(3).seedSelection(SeedSelection.ALL_KMERS)
                .queryPadding(2).graphPadding(1).maxSeedOccurrences(1)
                .minimumChainedSeeds(2)
                .fallback(SeededBandFallback.DENSE_IF_NO_CHAIN).build();
        SeededPOAAlignmentResult result = new SeededBandedPOAAligner(graph, scores)
                .align("AAAAAAAAAA", options);
        assertTrue(result.stats().repetitiveSeedsFiltered() > 0);
        assertTrue(result.stats().denseFallbackUsed());
    }

    @Test
    public void fullWidthSeedCorridorMatchesExistingDensePathExactly() {
        POGraphGlobal graph = new POGraphGlobal(Arrays.asList(
                "ACG-TACG", "ATGGTACG", "ACGGTTCG"), scores);
        String query = "ATGGTACG";
        List<POGNode> dense = graph.alignToGraph(query);
        SeededPOAAlignmentResult seeded = new SeededBandedPOAAligner(graph, scores)
                .align(query, options(3, 100, SeededBandFallback.NEVER));
        assertEquals(dense, seeded.path());
        assertEquals(seeded.stats().denseCellCount(), seeded.stats().evaluatedCells());
    }

    @Test
    public void localSeedCorridorFindsAnInternalQueryAlignment() {
        POGraphLocal graph = new POGraphLocal("ACGTTGCA", scores);
        SeededPOAAlignmentResult result = new SeededBandedPOAAligner(graph, scores)
                .align("TTACGTTGCAAA", options(3, 2, SeededBandFallback.NEVER));
        assertEquals("ACGTTGCA", pathBases(result.path()));
        assertEquals(2, result.queryInterval().start());
        assertEquals(10, result.queryInterval().end());
        assertTrue(result.stats().evaluatedCells() < result.stats().denseCellCount());
    }

    @Test
    public void fullWidthLocalCorridorMatchesDensePathScoreAndInterval() {
        POGraphLocal graph = new POGraphLocal(Arrays.asList(
                "ACG-TACG", "ATGGTACG", "ACGGTTCG"), scores);
        String query = "GGATGGTACGTT";
        List<POGNode> dense = graph.alignToGraph(query);
        jseqa.util.Interval denseInterval = graph.queryInterval();
        int denseScore = graph.alignToGraph(query, 1).get(0).score();
        SeededPOAAlignmentResult seeded = new SeededBandedPOAAligner(graph, scores)
                .align(query, options(3, 100, SeededBandFallback.NEVER));
        assertEquals(dense, seeded.path());
        assertEquals(denseScore, seeded.score());
        assertEquals(denseInterval.toString(), seeded.queryInterval().toString());
    }

    @Test
    public void reusableWorkspaceIsResetBetweenDifferentQueries() {
        POGraphLocal graph = new POGraphLocal("ACGTTGCA", scores);
        SeededBandedPOAAligner reused = new SeededBandedPOAAligner(graph, scores);
        SeededBandOptions options = options(3, 2, SeededBandFallback.NEVER);
        String shortQuery = "TTACGTTGCAAA";
        String longQuery = "GGGGACGTTGCATTTT";

        reused.align(shortQuery, options);
        SeededPOAAlignmentResult reusedLong = reused.align(longQuery, options);
        SeededPOAAlignmentResult reusedShort = reused.align(shortQuery, options);
        SeededPOAAlignmentResult freshLong =
                new SeededBandedPOAAligner(graph, scores).align(longQuery, options);
        SeededPOAAlignmentResult freshShort =
                new SeededBandedPOAAligner(graph, scores).align(shortQuery, options);

        assertEquals(freshLong.path(), reusedLong.path());
        assertEquals(freshLong.score(), reusedLong.score());
        assertEquals(freshLong.queryInterval().toString(),
                reusedLong.queryInterval().toString());
        assertEquals(freshShort.path(), reusedShort.path());
        assertEquals(freshShort.score(), reusedShort.score());
        assertEquals(freshShort.queryInterval().toString(),
                reusedShort.queryInterval().toString());
    }

    @Test
    public void localNoSeedFallbackMatchesExistingDenseAlignment() {
        POGraphLocal graph = new POGraphLocal("AAAAAAAA", scores);
        String query = "TTCCCCCCGG";
        List<POGNode> dense = graph.alignToGraph(query);
        jseqa.util.Interval interval = graph.queryInterval();
        SeededPOAAlignmentResult result = new SeededBandedPOAAligner(graph, scores)
                .align(query, options(4, 2, SeededBandFallback.DENSE_IF_NO_CHAIN));
        assertEquals(dense, result.path());
        assertEquals(interval.toString(), result.queryInterval().toString());
        assertTrue(result.stats().denseFallbackUsed());
    }

    @Test(expected = UnsupportedOperationException.class)
    public void fittingGraphModeIsRejectedExplicitly() {
        new SeededBandedPOAAligner(new POGraphFitting("ACGT", scores), scores);
    }

    private static SeededBandOptions options(int k, int padding,
            SeededBandFallback fallback) {
        return SeededBandOptions.builder()
                .kmerLength(k)
                .seedSelection(SeedSelection.ALL_KMERS)
                .queryPadding(padding)
                .graphPadding(1)
                .minimumChainedSeeds(2)
                .fallback(fallback)
                .build();
    }

    private static String pathBases(List<POGNode> path) {
        return path.stream().map(node -> node == null ? "-" : "" + node.getCH())
                .collect(Collectors.joining());
    }

    private static final class CountingPOGraphGlobal extends POGraphGlobal {
        int graphReads;

        CountingPOGraphGlobal(String sequence, CustomDnaSubstitutionMatrix scores) {
            super(sequence, scores);
        }

        @Override
        public org.jgrapht.graph.DefaultDirectedGraph<POGNode,POGEdge> getGraph() {
            graphReads++;
            return super.getGraph();
        }
    }
}
