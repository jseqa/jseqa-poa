package jseqa.poa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import jseqa.alignment.substitution.CustomDnaSubstitutionMatrix;
import org.junit.Test;

public class PositionedPOGraphTest {
    private static final PositionedPOGNodeFactory POSITIONED = new PositionedPOGNodeFactory();

    @Test
    public void initialSequenceRecordsEveryOriginalPosition() {
        POGraphFitting graph = new POGraphFitting("ACGT", POSITIONED);
        List<PositionedPOGNode> nodes = orderedNodes(graph);
        for (int position = 0; position < nodes.size(); position++) {
            assertEquals(position, nodes.get(position).positionForRead(0).getAsInt());
        }
    }

    @Test
    public void exactMatchAddsPositionsToExistingNodes() {
        POGraphFitting graph = new POGraphFitting("ACGT", POSITIONED);
        graph.add("ACGT");
        List<PositionedPOGNode> nodes = orderedNodes(graph);
        for (int position = 0; position < nodes.size(); position++) {
            assertEquals(position, nodes.get(position).positionForRead(1).getAsInt());
            assertEquals(2, nodes.get(position).getNumReads());
        }
    }

    @Test
    public void mismatchNodeRecordsItsQueryPosition() {
        POGraphFitting graph = new POGraphFitting("ACGT", mismatchPreferred(), POSITIONED);
        graph.add("ATGT");
        PositionedPOGNode mismatch = findObservation(graph, 1, 1, 'T');
        assertNotNull(mismatch);
        assertEquals(1, mismatch.positionForRead(1).getAsInt());
        assertFalse(mismatch.positionForRead(0).isPresent());
    }

    @Test
    public void insertionNodeRecordsItsQueryPosition() {
        POGraphFitting graph = new POGraphFitting("ACGT", gapPreferred(), POSITIONED);
        graph.add("ACAGT");
        PositionedPOGNode insertion = findObservation(graph, 1, 2, 'A');
        assertNotNull(insertion);
        assertFalse(insertion.positionForRead(0).isPresent());
    }

    @Test
    public void deletionAddsNoObservationForDeletedReferenceNode() {
        POGraphFitting graph = new POGraphFitting("ACGT", POSITIONED);
        graph.add("AGT");
        PositionedPOGNode deletedC = orderedNodes(graph).stream()
                .filter(node -> node.getCH() == 'C').findFirst().get();
        assertFalse(deletedC.positionForRead(1).isPresent());
        assertNotNull(findObservation(graph, 1, 0, 'A'));
        assertNotNull(findObservation(graph, 1, 1, 'G'));
        assertNotNull(findObservation(graph, 1, 2, 'T'));
    }

    @Test
    public void localAlignmentUsesOriginalQueryOffsets() {
        POGraphLocal graph = new POGraphLocal("ACGT", mismatchPreferred(), POSITIONED);
        graph.add("TTACGTGG");
        assertNotNull(findObservation(graph, 1, 2, 'A'));
        assertNotNull(findObservation(graph, 1, 3, 'C'));
        assertNotNull(findObservation(graph, 1, 4, 'G'));
        assertNotNull(findObservation(graph, 1, 5, 'T'));
    }

    @Test
    public void graphCopyPreservesIndependentPositionMaps() {
        POGraphFitting original = new POGraphFitting("ACGT", POSITIONED);
        POGraphFitting copy = new POGraphFitting(original);
        copy.add("ACGT");
        for (PositionedPOGNode node : orderedNodes(original)) {
            assertEquals(1, node.readPositions().size());
        }
        for (PositionedPOGNode node : orderedNodes(copy)) {
            assertEquals(2, node.readPositions().size());
        }
    }

    @Test
    public void everyAlignmentModeAcceptsPositionedFactory() {
        List<AbstractPartialOrderGraph> graphs = new ArrayList<AbstractPartialOrderGraph>();
        graphs.add(new POGraphGlobal("AC", POSITIONED));
        graphs.add(new POGraphLocal("AC", POSITIONED));
        graphs.add(new POGraphFitting("AC", POSITIONED));
        graphs.add(new POGraphOverlap("AC", POSITIONED));
        for (AbstractPartialOrderGraph graph : graphs) {
            assertTrue(graph.getGraph().vertexSet().stream()
                    .allMatch(node -> node instanceof PositionedPOGNode));
        }
    }

    private static PositionedPOGNode findObservation(AbstractPartialOrderGraph graph,
            int readIndex, int position, char residue) {
        return graph.getGraph().vertexSet().stream()
                .filter(node -> node instanceof PositionedPOGNode)
                .map(node -> (PositionedPOGNode) node)
                .filter(node -> node.getCH() == residue)
                .filter(node -> node.positionForRead(readIndex).isPresent())
                .filter(node -> node.positionForRead(readIndex).getAsInt() == position)
                .findFirst().orElse(null);
    }

    private static List<PositionedPOGNode> orderedNodes(AbstractPartialOrderGraph graph) {
        return graph.getGraph().vertexSet().stream()
                .map(node -> (PositionedPOGNode) node)
                .sorted(Comparator.comparingInt(node -> node.positionForRead(0).orElse(Integer.MAX_VALUE)))
                .collect(Collectors.toList());
    }

    private static CustomDnaSubstitutionMatrix mismatchPreferred() {
        return new CustomDnaSubstitutionMatrix(5, -1, -10, -1);
    }

    private static CustomDnaSubstitutionMatrix gapPreferred() {
        return new CustomDnaSubstitutionMatrix(5, -10, -2, -1);
    }
}
