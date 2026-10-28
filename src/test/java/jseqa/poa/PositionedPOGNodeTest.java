package jseqa.poa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class PositionedPOGNodeTest {

    @Test
    public void constructorRecordsReadAndPosition() {
        PositionedPOGNode node = new PositionedPOGNode('A', 3, 17);
        assertEquals(1, node.getNumReads());
        assertEquals(Integer.valueOf(17), node.readPositions().get(3));
        assertEquals(17, node.positionForRead(3).getAsInt());
    }

    @Test
    public void sameReadAndPositionIsIdempotent() {
        PositionedPOGNode node = new PositionedPOGNode('C', 2, 4);
        node.addReadPosition(2, 4);
        assertEquals(1, node.getNumReads());
        assertEquals(1, node.readPositions().size());
    }

    @Test(expected = IllegalStateException.class)
    public void sameReadCannotHaveTwoPositionsInOneNode() {
        PositionedPOGNode node = new PositionedPOGNode('G', 2, 4);
        node.addReadPosition(2, 5);
    }

    @Test
    public void missingReadReturnsEmptyPosition() {
        PositionedPOGNode node = new PositionedPOGNode('T', 0, 0);
        assertFalse(node.positionForRead(9).isPresent());
        assertTrue(node.positionForRead(0).isPresent());
    }

    @Test(expected = UnsupportedOperationException.class)
    public void readIndexCannotBeAddedWithoutPosition() {
        new PositionedPOGNode('A', 0, 0).addIDX(1);
    }
}
