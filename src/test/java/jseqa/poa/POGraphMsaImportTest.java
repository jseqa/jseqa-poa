package jseqa.poa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.junit.Test;

public class POGraphMsaImportTest {
    private static final List<String> MSA = Arrays.asList(
            "AC-GT",
            "ACTGT",
            "AC-GC");

    @Test
    public void fittingGraphPreservesAlignmentExactly() {
        POGraphFitting graph = new POGraphFitting(MSA);
        assertEquals(MSA, graph.getAlignment());
        assertEquals(6, graph.numNodes());
    }

    @Test
    public void importerMergesEqualResiduesAndAlignsAlternatives() {
        POGraphFitting graph = new POGraphFitting(MSA);
        POGNode sharedA = graph.getGraph().vertexSet().stream()
                .filter(node -> node.getCH() == 'A' && node.getNumReads() == 3)
                .findFirst().orElse(null);
        POGNode terminalT = graph.getGraph().vertexSet().stream()
                .filter(node -> node.getCH() == 'T' && node.getIDX().contains(0))
                .findFirst().orElse(null);
        assertNotNull(sharedA);
        assertNotNull(terminalT);
        assertTrue(terminalT.getIDX().contains(1));
        assertEquals(5, graph.getAlignment().get(0).length());
    }

    @Test
    public void edgesConnectReadsAcrossGapRuns() {
        POGraphFitting graph = new POGraphFitting(MSA);
        POGNode sharedC = findNode(graph, 'C', 0, 1);
        POGNode sharedG = findNode(graph, 'G', 0, 2);
        POGEdge skipGap = graph.getGraph().getEdge(sharedC, sharedG);
        assertNotNull(skipGap);
        assertTrue(graph.getEdgeIdx(skipGap).contains(0));
        assertTrue(graph.getEdgeIdx(skipGap).contains(2));
    }

    @Test
    public void positionedImportUsesUngappedReadCoordinates() {
        POGraphFitting graph = new POGraphFitting(MSA, new PositionedPOGNodeFactory());
        assertPosition(graph, 'T', 1, 2);
        assertPosition(graph, 'G', 0, 2);
        assertPosition(graph, 'G', 1, 3);
        assertPosition(graph, 'C', 2, 3);
    }

    @Test
    public void importedGraphCanAcceptAdditionalReads() {
        POGraphFitting graph = new POGraphFitting(MSA, new PositionedPOGNodeFactory());
        graph.add("ACGTC");
        assertEquals(4, graph.getAlignment().size());
        assertTrue(graph.getGraph().vertexSet().stream()
                .filter(node -> node instanceof PositionedPOGNode)
                .map(node -> (PositionedPOGNode) node)
                .anyMatch(node -> node.positionForRead(3).isPresent()));
    }

    @Test
    public void everyAlignmentModeCanBeBuiltFromMsa() {
        List<AbstractPartialOrderGraph> graphs = new ArrayList<AbstractPartialOrderGraph>();
        graphs.add(new POGraphGlobal(MSA));
        graphs.add(new POGraphLocal(MSA));
        graphs.add(new POGraphFitting(MSA));
        graphs.add(new POGraphOverlap(MSA));
        for (AbstractPartialOrderGraph graph : graphs) {
            assertEquals(MSA, graph.getAlignment());
        }
    }

    @Test(expected = IllegalArgumentException.class)
    public void rejectsRowsWithDifferentWidths() {
        new POGraphFitting(Arrays.asList("AC-G", "ACG"));
    }

    @Test(expected = IllegalArgumentException.class)
    public void rejectsAllGapColumns() {
        new POGraphFitting(Arrays.asList("A-C", "T-G"));
    }

    @Test(expected = IllegalArgumentException.class)
    public void rejectsAllGapRows() {
        new POGraphFitting(Arrays.asList("ACG", "---"));
    }

    @Test(expected = IllegalArgumentException.class)
    public void rejectsEmptyAlignment() {
        new POGraphFitting(Collections.<String>emptyList());
    }

    private static POGNode findNode(AbstractPartialOrderGraph graph, char residue,
            int readIndex, int readPosition) {
        return graph.getGraph().vertexSet().stream()
                .filter(node -> node.getCH() == residue)
                .filter(node -> node instanceof PositionedPOGNode
                        ? ((PositionedPOGNode) node).positionForRead(readIndex).orElse(-1) == readPosition
                        : node.getIDX().contains(readIndex))
                .findFirst().orElse(null);
    }

    private static void assertPosition(AbstractPartialOrderGraph graph, char residue,
            int readIndex, int expectedPosition) {
        PositionedPOGNode node = graph.getGraph().vertexSet().stream()
                .filter(candidate -> candidate instanceof PositionedPOGNode)
                .map(candidate -> (PositionedPOGNode) candidate)
                .filter(candidate -> candidate.getCH() == residue)
                .filter(candidate -> candidate.positionForRead(readIndex).orElse(-1) == expectedPosition)
                .findFirst().orElse(null);
        assertNotNull(node);
    }
}
