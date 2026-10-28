package jseqa.poa;

/** Creates {@link PositionedPOGNode} instances for a partial-order graph. */
public final class PositionedPOGNodeFactory implements POGNodeFactory {

    @Override
    public POGNode create(char residue, int readIndex, int readPosition) {
        return new PositionedPOGNode(residue, readIndex, readPosition);
    }

    @Override
    public POGNode copy(POGNode node) {
        if (!(node instanceof PositionedPOGNode)) {
            throw new IllegalArgumentException("cannot copy a node without position metadata");
        }
        return new PositionedPOGNode((PositionedPOGNode) node);
    }

    @Override
    public void addReadObservation(POGNode node, int readIndex, int readPosition) {
        if (!(node instanceof PositionedPOGNode)) {
            throw new IllegalArgumentException("node was not created by this factory");
        }
        ((PositionedPOGNode) node).addReadPosition(readIndex, readPosition);
    }
}
