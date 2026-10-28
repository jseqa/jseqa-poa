package jseqa.poa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;
import jseqa.alignment.substitution.CustomDnaSubstitutionMatrix;
import org.junit.Test;

public class POGraphSourceRecurrenceTest {

    @Test
    public void sourceMismatchCreatesAlternativeNode() {
        POGraphFitting graph = new POGraphFitting("A", mismatchPreferredMatrix());
        graph.add("T");
        assertEquals(2, graph.numNodes());
        assertEquals("A", graph.getAlignment().get(0));
        assertEquals("T", graph.getAlignment().get(1));
    }

    @Test
    public void globalAlignmentRetainsLeadingQueryInsertion() {
        List<POGNode> path = new POGraphGlobal("A").alignToGraph("TA");
        assertEquals(2, path.size());
        assertNull(path.get(0));
        assertEquals('A', path.get(1).getCH());
    }

    @Test
    public void secondGraphSourceBacktracksToOrigin() {
        POGraphFitting graph = new POGraphFitting("AC", mismatchPreferredMatrix());
        graph.add("TC");
        List<POGNode> path = graph.alignToGraph("TC");
        assertEquals(2, path.size());
        assertEquals('T', path.get(0).getCH());
        assertEquals('C', path.get(1).getCH());
    }

    private static CustomDnaSubstitutionMatrix mismatchPreferredMatrix() {
        return new CustomDnaSubstitutionMatrix(5, -1, -10, -1);
    }
}
