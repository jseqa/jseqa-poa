package jseqa.poa;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class POGraphMatrixSizingTest {

    @Test
    public void acceptsQueryAtFormerColumnLimit() {
        String sequence = repeat('A', 500);
        assertEquals(500, new POGraphFitting(sequence).alignToGraph(sequence).size());
    }

    @Test
    public void acceptsQueryBeyondFormerColumnLimit() {
        String sequence = repeat('C', 501);
        assertEquals(501, new POGraphGlobal(sequence).alignToGraph(sequence).size());
    }

    @Test
    public void acceptsGraphAtFormerRowLimit() {
        POGraphLocal graph = new POGraphLocal(repeat('A', 2000));
        assertFalse(graph.alignToGraph("A").isEmpty());
    }

    private static String repeat(char value, int count) {
        StringBuilder result = new StringBuilder(count);
        for (int i = 0; i < count; i++) {
            result.append(value);
        }
        return result.toString();
    }
}
