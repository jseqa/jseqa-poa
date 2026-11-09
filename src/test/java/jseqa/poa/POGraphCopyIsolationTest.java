package jseqa.poa;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;
import org.junit.Test;

public class POGraphCopyIsolationTest {

    @Test
    public void addingToCopyDoesNotChangeOriginalNodeCount() {
        POGraphFitting original = new POGraphFitting("AC");
        POGraphFitting copy = new POGraphFitting(original);
        copy.add("AT");
        assertEquals(2, original.numNodes());
        assertEquals(3, copy.numNodes());
    }

    @Test
    public void addingToCopyDoesNotChangeOriginalReadMembership() {
        POGraphOverlap original = new POGraphOverlap("AC");
        List<Integer> before = new ArrayList<Integer>();
        for (POGNode node : original.getGraph().vertexSet()) {
            before.add(node.getNumReads());
        }
        new POGraphOverlap(original).add("AC");
        List<Integer> after = new ArrayList<Integer>();
        for (POGNode node : original.getGraph().vertexSet()) {
            after.add(node.getNumReads());
        }
        assertEquals(before, after);
    }

    @Test
    public void copiedGlobalGraphInitializesEndpointState() {
        POGraphGlobal copy = new POGraphGlobal(new POGraphFitting("AC"));
        assertEquals(2, copy.alignToGraph("AC").size());
        copy.add("AC");
        assertEquals(2, copy.getAlignment().size());
    }
}
