package jseqa.poa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.List;
import jseqa.alignment.substitution.CustomDnaSubstitutionMatrix;
import org.junit.Test;

public class GlobalProfilePOAMergerTest {
    private final CustomDnaSubstitutionMatrix scores =
            new CustomDnaSubstitutionMatrix(5, -4, -4, -1);
    private final GlobalProfilePOAMerger merger = new GlobalProfilePOAMerger(scores);

    @Test
    public void mergesSequenceWithSequenceAndPreservesMismatchBranch() {
        PartialOrderAlignment merged = merger.merge(
                new POGraphGlobal("ACGT", scores),
                new POGraphGlobal("AGGT", scores));
        assertReconstructs(Arrays.asList("ACGT", "AGGT"), merged.getAlignment());
        assertTrue(merged.numNodes() >= 5);
    }

    @Test
    public void mergesSingleSequenceWithMultiSequenceProfile() {
        PartialOrderAlignment left = new POGraphGlobal("ACGTT", scores);
        PartialOrderAlignment right = new POGraphGlobal(
                Arrays.asList("AC-GT", "ATCGT"), scores);
        PartialOrderAlignment merged = merger.merge(left, right);
        assertReconstructs(Arrays.asList("ACGTT", "ACGT", "ATCGT"),
                merged.getAlignment());
    }

    @Test
    public void mergesTwoMultiSequenceProfilesWithIndels() {
        PartialOrderAlignment left = new POGraphGlobal(
                Arrays.asList("-ACGT-", "AACGTT"), scores);
        PartialOrderAlignment right = new POGraphGlobal(
                Arrays.asList("AC-GTT", "ACCGTT"), scores);
        PartialOrderAlignment merged = merger.merge(left, right);
        assertReconstructs(Arrays.asList("ACGT", "AACGTT", "ACGTT", "ACCGTT"),
                merged.getAlignment());
    }

    @Test
    public void duplicateSequencesRemainDistinctColors() {
        PartialOrderAlignment merged = merger.merge(
                new POGraphGlobal(Arrays.asList("ACGT", "ACGT"), scores),
                new POGraphGlobal("ACGT", scores));
        assertEquals(3, merged.getAlignment().size());
        assertReconstructs(Arrays.asList("ACGT", "ACGT", "ACGT"),
                merged.getAlignment());
        for (POGNode node : merged.getGraph().vertexSet()) {
            assertEquals(3, node.getNumReads());
        }
    }

    private static void assertReconstructs(List<String> expected, List<String> rows) {
        assertEquals(expected.size(), rows.size());
        int width = rows.get(0).length();
        for (int i = 0; i < rows.size(); i++) {
            assertEquals(width, rows.get(i).length());
            assertEquals(expected.get(i), rows.get(i).replace("-", ""));
        }
    }
}
