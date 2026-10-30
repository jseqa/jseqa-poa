package jseqa.poa;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jseqa.alignment.substitution.CustomDnaSubstitutionMatrix;
import jseqa.alignment.substitution.SubstitutionMatrix;
import org.junit.Test;

public class POGraphProgressiveMsaPreservationTest {
    private static final List<String> MINIMAL = Arrays.asList(
            "TGCGAGA",
            "CAACAGA");

    @Test
    public void globalProgressiveMsaPreservesMinimalReproducer() {
        assertPreserved(new POGraphGlobal(scoring()), MINIMAL);
    }

    @Test
    public void everyAlignmentModePreservesMinimalReproducer() {
        for (AbstractPartialOrderGraph graph : allModes()) {
            assertPreserved(graph, MINIMAL);
        }
    }

    @Test
    public void repeatedMotifsRemainOrderedForSeveralInputOrders() {
        List<List<String>> orders = Arrays.asList(
                Arrays.asList("TGCGAGA", "CAACAGA", "TGACAGA"),
                Arrays.asList("CAACAGA", "TGCGAGA", "TGACAGA"),
                Arrays.asList("TGACAGA", "CAACAGA", "TGCGAGA"));

        for (List<String> sequences : orders) {
            for (AbstractPartialOrderGraph graph : allModes()) {
                assertPreserved(graph, sequences);
            }
        }
    }

    @Test
    public void progressiveMsaIsDeterministicAcrossRepeatedBuilds() {
        List<String> expected = null;
        for (int run = 0; run < 10; run++) {
            POGraphGlobal graph = new POGraphGlobal(scoring());
            graph.addAll(MINIMAL);
            List<String> alignment = graph.getAlignment();
            assertUngappedEquals("deterministic global build", MINIMAL, alignment);
            if (expected == null) {
                expected = alignment;
            } else {
                assertEquals(expected, alignment);
            }
        }
    }

    @Test(expected = IllegalStateException.class)
    public void duplicateReadMembershipWithinAColumnIsRejected() {
        POGraphGlobal graph = new POGraphGlobal("AC");
        POGNode first = null;
        POGNode second = null;
        for (POGNode node : graph.g.vertexSet()) {
            if (node.getCH() == 'A') {
                first = node;
            } else if (node.getCH() == 'C') {
                second = node;
            }
        }
        graph.alignedG.addEdge(first, second);
        graph.getAlignment();
    }

    @Test(expected = IllegalStateException.class)
    public void cyclicAlignmentColumnQuotientIsRejected() {
        POGraphGlobal graph = new POGraphGlobal();
        POGNode a = new POGNode('A', 0);
        POGNode b = new POGNode('B', 0);
        POGNode c = new POGNode('C', 1);
        POGNode d = new POGNode('D', 1);
        for (POGNode node : Arrays.asList(a, b, c, d)) {
            graph.g.addVertex(node);
            graph.alignedG.addVertex(node);
        }
        graph.g.addEdge(a, b);
        graph.g.addEdge(c, d);
        graph.alignedG.addEdge(a, d);
        graph.alignedG.addEdge(b, c);
        graph.readIDX = 2;

        graph.getAlignment();
    }

    private static List<AbstractPartialOrderGraph> allModes() {
        SubstitutionMatrix score = scoring();
        List<AbstractPartialOrderGraph> graphs = new ArrayList<AbstractPartialOrderGraph>();
        graphs.add(new POGraphGlobal(score));
        graphs.add(new POGraphFitting(score));
        graphs.add(new POGraphOverlap(score));
        graphs.add(new POGraphLocal(score));
        return graphs;
    }

    private static SubstitutionMatrix scoring() {
        return new CustomDnaSubstitutionMatrix(5, -4, -4, -1);
    }

    private static void assertPreserved(AbstractPartialOrderGraph graph,
            List<String> sequences) {
        graph.addAll(sequences);
        assertUngappedEquals(graph.getClass().getSimpleName(), sequences, graph.getAlignment());
    }

    private static void assertUngappedEquals(String message, List<String> expected,
            List<String> alignment) {
        List<String> ungapped = new ArrayList<String>(alignment.size());
        for (String row : alignment) {
            ungapped.add(row.replace("-", ""));
        }
        assertEquals(message, expected, ungapped);
    }
}
