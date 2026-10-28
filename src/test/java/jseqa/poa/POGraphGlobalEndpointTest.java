package jseqa.poa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import org.junit.Test;

public class POGraphGlobalEndpointTest {

    @Test
    public void singleNodeReadIsBothStartAndEnd() {
        POGraphGlobal graph = new POGraphGlobal("A");
        graph.add("A");
        assertEquals(1, graph.numNodes());
        assertEquals(2, graph.getGraph().vertexSet().iterator().next().getNumReads());
    }

    @Test
    public void backtrackingIncludesLastNode() {
        List<POGNode> path = new POGraphGlobal("AC").alignToGraph("AC");
        assertEquals(2, path.size());
        assertEquals('C', path.get(1).getCH());
    }

    @Test
    public void terminalMismatchCreatesSupportedAlternativeEnd() {
        POGraphGlobal graph = new POGraphGlobal("AC");
        graph.add("AT");
        POGNode terminalT = graph.getGraph().vertexSet().stream()
                .filter(node -> node.getCH() == 'T').findFirst().orElse(null);
        assertNotNull(terminalT);
        assertTrue(terminalT.getIDX().contains(1));
        assertEquals(1, graph.getGraph().inDegreeOf(terminalT));
    }
}
