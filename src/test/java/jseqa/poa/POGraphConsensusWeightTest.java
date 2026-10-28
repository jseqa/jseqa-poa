package jseqa.poa;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class POGraphConsensusWeightTest {

    @Test
    public void consensusCountsSupportOnFirstNode() {
        assertSupportedSourceWins(5, 1, 2, "HE");
    }

    @Test
    public void consensusCountsSupportAcrossThreeWayDifference() {
        assertSupportedSourceWins(8, 2, 4, "HE");
    }

    @Test
    public void consensusCountsSourceWhenMarginIsOne() {
        assertSupportedSourceWins(4, 1, 3, "HE");
    }

    private static void assertSupportedSourceWins(int highSupport, int lowSupport,
            int intermediateSupport, String expected) {
        POGraphFitting graph = new POGraphFitting();
        POGNode high = node('H', highSupport);
        POGNode low = node('L', lowSupport);
        POGNode intermediate = node('I', intermediateSupport);
        POGNode end = node('E', 1);
        graph.getGraph().addVertex(high);
        graph.getGraph().addVertex(low);
        graph.getGraph().addVertex(intermediate);
        graph.getGraph().addVertex(end);
        graph.getGraph().addEdge(high, end);
        graph.getGraph().addEdge(low, intermediate);
        graph.getGraph().addEdge(intermediate, end);
        assertEquals(expected, graph.consensus());
    }

    private static POGNode node(char residue, int support) {
        POGNode node = new POGNode(residue, 0);
        for (int i = 1; i < support; i++) {
            node.addIDX(i);
        }
        return node;
    }
}
