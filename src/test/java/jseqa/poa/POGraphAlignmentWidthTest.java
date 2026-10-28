package jseqa.poa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.Arrays;
import java.util.List;
import org.junit.Test;

public class POGraphAlignmentWidthTest {

    @Test
    public void fittingAlignmentHasNoPaddingColumns() {
        assertExactWidth(new POGraphFitting(), "AC", "AT");
    }

    @Test
    public void overlapAlignmentHasNoPaddingColumns() {
        assertExactWidth(new POGraphOverlap(), "ACG", "AG");
    }

    @Test
    public void globalAlignmentHasNoPaddingColumns() {
        assertExactWidth(new POGraphGlobal(), "AC", "TAC");
    }

    private static void assertExactWidth(AbstractPartialOrderGraph graph, String... sequences) {
        graph.addAll(Arrays.asList(sequences));
        int columns = 0;
        for (POGNode ignored : graph) {
            columns++;
        }
        List<String> alignment = graph.getAlignment();
        for (String row : alignment) {
            assertEquals(columns, row.length());
            assertFalse(row.endsWith("----------"));
        }
    }
}
