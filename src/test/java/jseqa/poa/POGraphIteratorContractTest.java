package jseqa.poa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Iterator;
import jseqa.alignment.substitution.CustomDnaSubstitutionMatrix;
import org.junit.Test;

public class POGraphIteratorContractTest {

    @Test
    public void iteratorHandlesAlternativeInLastColumn() {
        POGraphFitting graph = new POGraphFitting("AC", mismatchPreferredMatrix());
        graph.add("AT");
        assertEquals(2, count(graph.iterator()));
    }

    @Test
    public void iteratorHandlesOnlyAlternativeColumn() {
        POGraphFitting graph = new POGraphFitting("A", mismatchPreferredMatrix());
        graph.add("T");
        assertEquals(1, count(graph.iterator()));
    }

    @Test
    public void hasNextIsFalseAfterLastUnvisitedColumn() {
        POGraphOverlap graph = new POGraphOverlap("AC", mismatchPreferredMatrix());
        graph.add("AT");
        Iterator<POGNode> iterator = graph.iterator();
        assertTrue(iterator.hasNext());
        iterator.next();
        assertTrue(iterator.hasNext());
        iterator.next();
        assertFalse(iterator.hasNext());
    }

    private static int count(Iterator<POGNode> iterator) {
        int count = 0;
        while (iterator.hasNext()) {
            iterator.next();
            count++;
        }
        return count;
    }

    private static CustomDnaSubstitutionMatrix mismatchPreferredMatrix() {
        return new CustomDnaSubstitutionMatrix(5, -1, -10, -1);
    }
}
