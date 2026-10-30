package jseqa.poa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import jseqa.alignment.AlignmentMatrixScore;
import jseqa.poa.POGPtr.PointerType;
import org.junit.Test;

public class POGraphPerformanceCompatibilityTest {

    @Test
    public void reusesMatricesForLongThenShortQuery() {
        POGraphLocal graph = new POGraphLocal("ACGTACGTACGT");
        graph.alignToGraph("TTACGTACGTACGTGG");
        AlignmentMatrixScore scores = graph.alignScoreMat;
        AlignPOAMatrixPtr coordinatePointers = graph.alignPtrMat;
        AlignPOGPtrMatrix operationPointers = graph.ptrMat;

        assertEquals("ACGT", pathString(graph.alignToGraph("ACGT")));
        assertSame(scores, graph.alignScoreMat);
        assertSame(coordinatePointers, graph.alignPtrMat);
        assertSame(operationPointers, graph.ptrMat);
        assertEquals(5, graph.alignScoreMat.ncols());
    }

    @Test
    public void reusesMatricesForShortThenLongQuery() {
        POGraphFitting graph = new POGraphFitting("ACGTACGT");
        graph.alignToGraph("ACGT");
        AlignmentMatrixScore scores = graph.alignScoreMat;

        assertEquals("ACGTACGT", pathString(graph.alignToGraph("TTACGTACGTAA")));
        assertSame(scores, graph.alignScoreMat);
        assertEquals(13, graph.alignScoreMat.ncols());
    }

    @Test
    public void repeatedIdenticalQueriesReturnExactlySamePath() {
        POGraphGlobal graph = new POGraphGlobal("ACGTACGT");
        List<POGNode> first = graph.alignToGraph("ACGTACGT");
        List<POGNode> second = graph.alignToGraph("ACGTACGT");
        assertEquals(first, second);
    }

    @Test
    public void equalScoringPredecessorsKeepTheFirstCandidate() {
        TieProbeGraph graph = new TieProbeGraph();
        graph.prepareProbe();
        graph.alignScoreMat.set(1, 0, 7);
        graph.alignScoreMat.set(2, 0, 7);

        graph.probeCellWithPredecessors(3, 1, new int[] {2, 1});

        assertEquals(2, graph.ptrMat.previousRow(3, 1));
        assertEquals(0, graph.ptrMat.previousColumn(3, 1));
        assertEquals(PointerType.MATCH, graph.ptrMat.type(3, 1));
    }

	@Test
	public void cumulativeColorCoherenceResolvesEqualScoringPredecessors() {
		ColorTieProbeGraph graph = new ColorTieProbeGraph();
		graph.prepareProbe();
		graph.setColorCoherenceTieBreaking(true);
		int sharedRow = graph.rowFor('A');
		int narrowRow = graph.rowFor('T');
		int currentRow = graph.rowFor('C');
		graph.alignScoreMat.set(sharedRow, 0, 7);
		graph.alignScoreMat.set(narrowRow, 0, 7);

		graph.probeCellWithPredecessors(currentRow, 1,
				new int[] {narrowRow, sharedRow});

		assertEquals(sharedRow, graph.ptrMat.previousRow(currentRow, 1));
		assertEquals(PointerType.MATCH, graph.ptrMat.type(currentRow, 1));
	}

    @Test
    public void repeatedReadOnlyAlignmentsReuseTopologyCache() {
        POGraphLocal graph = new POGraphLocal("ACGTACGT");
        int before = graph.topologyBuildCount();
        graph.alignToGraph("ACGT");
        int afterFirst = graph.topologyBuildCount();
        graph.alignToGraph("CGTA");
        graph.alignToGraph("GTAC");
        assertTrue(afterFirst > before);
        assertEquals(afterFirst, graph.topologyBuildCount());
    }

    @Test
    public void graphMutationInvalidatesTopologyCache() {
        POGraphFitting graph = new POGraphFitting("ACGT");
        graph.alignToGraph("ACGT");
        int beforeMutation = graph.topologyBuildCount();
        graph.add("ACAT");
        assertEquals("ACAT", pathString(graph.alignToGraph("ACAT")));
		int afterMutation = graph.topologyBuildCount();

		assertTrue(afterMutation > beforeMutation);
        assertEquals(afterMutation, graph.topologyBuildCount());
    }

    @Test
    public void branchingGraphWithSeveralPredecessorsKeepsExactPaths() {
        List<String> alignment = Arrays.asList(
                "AC-GT",
                "AT-GT",
                "ACAGT",
                "ATAGT");
        POGraphFitting graph = new POGraphFitting(alignment);
        for (String row : alignment) {
            String query = row.replace("-", "");
            assertEquals(query, pathString(graph.alignToGraph(query)));
        }
    }

    @Test
    public void positionedNodesRetainCoordinatesWithReusedMatrices() {
        POGraphFitting graph = new POGraphFitting("ACGT",
                new PositionedPOGNodeFactory());
        graph.alignToGraph("ACGTACGT");
        graph.add("ACAT");
        PositionedPOGNode mismatch = graph.getGraph().vertexSet().stream()
                .filter(node -> node instanceof PositionedPOGNode)
                .map(node -> (PositionedPOGNode) node)
                .filter(node -> node.getCH() == 'A')
                .filter(node -> node.positionForRead(1).orElse(-1) == 2)
                .findFirst().orElse(null);
        assertTrue(mismatch != null);
    }

    private static String pathString(List<POGNode> path) {
        return path.stream()
                .filter(node -> node != null)
                .map(node -> Character.toString(node.getCH()))
                .collect(Collectors.joining());
    }

    private static final class TieProbeGraph extends POGraphGlobal {
        TieProbeGraph() {
            super("AAA");
        }

        void prepareProbe() {
            computeInternalDS();
            prepareMatrices(1);
        }

        void probeCellWithPredecessors(int row, int column, int[] predecessors) {
            fillCell(row, column, 'A', 'A', predecessors);
        }
    }

	private static final class ColorTieProbeGraph extends POGraphGlobal {
		ColorTieProbeGraph() {
			super(Arrays.asList("A-C", "A-C", "ATC"));
		}

		void prepareProbe() {
			computeInternalDS();
			prepareMatrices(1);
		}

		int rowFor(char base) {
			for (int index = 0; index < nodeList.size(); index++) {
				if (nodeList.get(index).getCH() == base) {
					return index + 1;
				}
			}
			throw new IllegalArgumentException("No node for " + base);
		}

		void probeCellWithPredecessors(int row, int column, int[] predecessors) {
			fillCell(row, column, 'C', 'C', predecessors);
		}
	}
}
