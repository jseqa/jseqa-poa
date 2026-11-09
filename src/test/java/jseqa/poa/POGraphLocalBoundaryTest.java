package jseqa.poa;

import static org.junit.Assert.assertEquals;

import java.util.List;
import org.junit.Test;

public class POGraphLocalBoundaryTest {

    @Test
    public void singleCharacterOptimumIncludesFinalRowAndColumn() {
        assertPath("A", "A", "A");
    }

    @Test
    public void optimumMayUseFinalQueryCharacter() {
        assertPath("AB", "ZA", "A");
    }

    @Test
    public void optimumMayUseFinalGraphNode() {
        assertPath("BA", "AZ", "A");
    }

    private static void assertPath(String graphSequence, String query, String expected) {
        List<POGNode> path = new POGraphLocal(graphSequence).alignToGraph(query);
        StringBuilder actual = new StringBuilder();
        for (POGNode node : path) {
            if (node != null) {
                actual.append(node.getCH());
            }
        }
        assertEquals(expected, actual.toString());
    }
}
