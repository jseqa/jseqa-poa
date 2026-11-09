package jseqa.poa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.List;

import org.junit.Test;

import jseqa.util.Interval;

/** Tests ranked sequence-to-graph path backtracking. */
public class POGraphKBestAlignmentTest {

    @Test(expected = IllegalArgumentException.class)
    public void rejectsNonPositivePathCount() {
        new POGraphLocal(Arrays.asList("ACG")).alignToGraph("ACG", 0);
    }

    @Test
    public void rankOneMatchesLegacyBacktracking() {
        POGraphLocal graph = new POGraphLocal(Arrays.asList("ACG", "ATG", "A-G"));

        List<POGNode> expectedPath = graph.alignToGraph("AAG");
        Interval expectedInterval = graph.queryInterval();
        POAAlignmentPath actual = graph.alignToGraph("AAG", 1).get(0);

        assertEquals(expectedPath, actual.path());
        assertEquals(expectedInterval.start(), actual.queryInterval().start());
        assertEquals(expectedInterval.end(), actual.queryInterval().end());
    }

    @Test
    public void bestAlignmentReturnsLegacyPathAndScoreInOnePass() {
        POGraphLocal graph = new POGraphLocal(Arrays.asList("ACGT", "ACCT"));
        String query = "TTACGTAA";
        List<POGNode> expectedPath = graph.alignToGraph(query);
        Interval expectedInterval = graph.queryInterval();
        int expectedScore = graph.alignToGraph(query, 1).get(0).score();

        POAAlignmentPath actual = graph.alignBestToGraph(query);

        assertEquals(expectedPath, actual.path());
        assertEquals(expectedInterval.start(), actual.queryInterval().start());
        assertEquals(expectedInterval.end(), actual.queryInterval().end());
        assertEquals(expectedScore, actual.score());
    }

    @Test
    public void returnsDistinctPathsInDescendingScoreOrder() {
        POGraphLocal graph = new POGraphLocal(Arrays.asList("ACG", "ATG", "A-G"));

        List<POAAlignmentPath> paths = graph.alignToGraph("AAG", 3);

        assertTrue(paths.size() > 1);
        assertTrue(paths.size() <= 3);
        for (int index = 1; index < paths.size(); index++) {
            assertTrue(paths.get(index - 1).score() >= paths.get(index).score());
            assertFalse(samePath(paths.get(index - 1), paths.get(index)));
        }
    }

    private static boolean samePath(POAAlignmentPath left, POAAlignmentPath right) {
        return left.queryInterval().start() == right.queryInterval().start()
                && left.queryInterval().end() == right.queryInterval().end()
                && left.path().equals(right.path());
    }
}
