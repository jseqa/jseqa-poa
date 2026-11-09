package jseqa.poa;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class POGraphEdgeBookkeepingTest {

    @Test
    public void terminalMatchDoesNotCreateNullEdgeEntry() {
        POGraphFitting graph = new POGraphFitting("AC");
        graph.add("AC");
        assertFalse(graph.edgeReadIdx.containsKey(null));
    }

    @Test
    public void terminalInsertionDoesNotCreateNullEdgeEntry() {
        POGraphGlobal graph = new POGraphGlobal("A");
        graph.add("AA");
        assertFalse(graph.edgeReadIdx.containsKey(null));
    }

    @Test
    public void everyAlternativeEdgeHasReadMetadata() {
        POGraphFitting graph = new POGraphFitting("ACA");
        graph.add("ATA");
        for (POGEdge edge : graph.getGraph().edgeSet()) {
            assertNotNull(graph.getEdgeIdx(edge));
        }
        assertTrue(graph.getGraph().edgeSet().stream()
                .anyMatch(edge -> graph.getEdgeIdx(edge).contains(1)));
    }
}
