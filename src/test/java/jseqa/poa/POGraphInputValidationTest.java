package jseqa.poa;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class POGraphInputValidationTest {

    @Test(expected = IllegalArgumentException.class)
    public void constructorRejectsEmptySequence() {
        new POGraphGlobal("");
    }

    @Test
    public void addRejectsEmptySequenceWithoutChangingGraph() {
        POGraphFitting graph = new POGraphFitting("AC");
        try {
            graph.add("");
        } catch (IllegalArgumentException expected) {
            assertEquals(2, graph.numNodes());
            assertEquals(1, graph.getAlignment().size());
            return;
        }
        throw new AssertionError("Expected empty sequence to be rejected");
    }

    @Test(expected = IllegalArgumentException.class)
    public void alignmentRejectsNullSequence() {
        new POGraphLocal("AC").alignToGraph(null);
    }
}
