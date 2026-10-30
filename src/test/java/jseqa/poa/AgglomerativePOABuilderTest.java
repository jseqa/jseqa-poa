package jseqa.poa;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import jseqa.alignment.substitution.CustomDnaSubstitutionMatrix;
import org.junit.Test;

public class AgglomerativePOABuilderTest {
    private final CustomDnaSubstitutionMatrix scores =
            new CustomDnaSubstitutionMatrix(5, -4, -4, -1);

    @Test
    public void pairwiseScoresAreSymmetricAndIdenticalSequencesAreMaximal() {
        double[][] matrix = new AgglomerativePOABuilder(scores)
                .pairwiseSimilarities(Arrays.asList("ACGT", "ACGT", "AGGT"));
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                assertEquals(matrix[i][j], matrix[j][i], 0.0);
            }
        }
        assertEquals(1.0, matrix[0][1], 0.0);
        assertTrue(matrix[0][1] > matrix[0][2]);
    }

    @Test
    public void parallelAndSingleThreadedPairwiseScoresAreIdentical() {
        List<String> input = Arrays.asList("ACGT", "AGGT", "ACGTT", "TTTT");
        double[][] single = new AgglomerativePOABuilder(scores).threads(1)
                .pairwiseSimilarities(input);
        double[][] parallel = new AgglomerativePOABuilder(scores).threads(3)
                .pairwiseSimilarities(input);
        for (int i = 0; i < single.length; i++) {
            assertArrayEquals(single[i], parallel[i], 0.0);
        }
    }

    @Test
    public void meanLengthNormalizationDividesRawScoreByMeanLength() {
        List<String> input = Arrays.asList("ACGT", "ACGTTT");
        double raw = new AgglomerativePOABuilder(
                scores, Linkage.UPGMA, ScoreNormalization.RAW)
                .pairwiseSimilarities(input)[0][1];
        double normalized = new AgglomerativePOABuilder(
                scores, Linkage.UPGMA, ScoreNormalization.MEAN_LENGTH)
                .pairwiseSimilarities(input)[0][1];
        assertEquals(raw / 5.0, normalized, 0.0);
    }

    @Test
    public void closestPairMergesFirstAndUpgmaUsesLeafPairAverage() {
        List<String> input = Arrays.asList("AAAA", "AAAT", "CCCC");
        AgglomerativePOABuilder builder = new AgglomerativePOABuilder(
                scores, Linkage.UPGMA, ScoreNormalization.SELF_SCORE);
        double[][] similarities = builder.pairwiseSimilarities(input);
        ProgressivePOAResult result = builder.build(input);
        GuideTreeMerge first = result.guideTree().merges().get(0);
        GuideTreeMerge second = result.guideTree().merges().get(1);
        assertEquals(0, first.leftCluster());
        assertEquals(1, first.rightCluster());
        assertEquals((similarities[0][2] + similarities[1][2]) / 2.0,
                second.similarity(), 0.0);
    }

    @Test
    public void tiesAndThreadCountsProduceDeterministicTreeAndAlignment() {
        List<String> input = Arrays.asList("ACGT", "ACGT", "ACGT", "ACGA");
        ProgressivePOAResult first = new AgglomerativePOABuilder(scores)
                .threads(1).build(input);
        ProgressivePOAResult second = new AgglomerativePOABuilder(scores)
                .threads(4).build(input);
        assertEquals(first.guideTree().toNewick(), second.guideTree().toNewick());
        assertEquals(first.alignment(), second.alignment());
        assertEquals(first.alignmentInputIndices(), second.alignmentInputIndices());
        assertEquals(0, first.guideTree().merges().get(0).leftCluster());
        assertEquals(1, first.guideTree().merges().get(0).rightCluster());
    }

    @Test
    public void everyInputOccurrenceIsMappedAndReconstructedExactlyOnce() {
        List<String> input = Arrays.asList(
                "TGCGAGA", "CAACAGA", "TGCGAGA", "TGACAGA", "AACAGAA");
        ProgressivePOAResult result = new AgglomerativePOABuilder(scores).build(input);
        assertEquals(input.size(), result.alignment().size());
        Set<Integer> mapped = new HashSet<Integer>(result.alignmentInputIndices());
        assertEquals(input.size(), mapped.size());
        for (int row = 0; row < result.alignment().size(); row++) {
            int inputIndex = result.alignmentInputIndices().get(row);
            assertEquals(input.get(inputIndex),
                    result.alignment().get(row).replace("-", ""));
        }
    }

    @Test
    public void singleSequenceSkipsClustering() {
        ProgressivePOAResult result = new AgglomerativePOABuilder(scores)
                .build(Arrays.asList("ACGT"));
        assertTrue(result.guideTree().root().isLeaf());
        assertTrue(result.guideTree().merges().isEmpty());
        assertEquals(Arrays.asList("ACGT"), result.alignment());
        assertEquals(Arrays.asList(0), result.alignmentInputIndices());
    }

    @Test
    public void ighv169TerminalInsertionKeepsSharedAgaSuffix() {
        String allele05 = "GTATTACTGTGCGAGA";
        String allele12 = "GTATTACTGTGCGAGAGA";
        ProgressivePOAResult result = new AgglomerativePOABuilder(scores)
                .build(Arrays.asList(allele05, allele12));
        List<String> rows = result.alignment();
        assertEquals(rows.get(0).length(), rows.get(1).length());
        assertTrue(rows.get(0).endsWith("AGA"));
        assertTrue(rows.get(1).endsWith("AGA"));
        assertEquals(allele05, rows.get(0).replace("-", ""));
        assertEquals(allele12, rows.get(1).replace("-", ""));
        assertEquals(2, count(rows.get(0), '-'));
    }

    @Test(expected = IllegalArgumentException.class)
    public void rejectsEmptyInput() {
        new AgglomerativePOABuilder(scores).build(Arrays.<String>asList());
    }

    private static int count(String value, char target) {
        int count = 0;
        for (int i = 0; i < value.length(); i++) {
            if (value.charAt(i) == target) count++;
        }
        return count;
    }
}
