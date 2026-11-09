package jseqa.poa;

/**
 * Creates graph nodes and records observations on them. Implementations can
 * attach metadata beyond the residue and read index retained by {@link POGNode}.
 */
public interface POGNodeFactory {

    POGNode create(char residue, int readIndex, int readPosition);

    POGNode copy(POGNode node);

    void addReadObservation(POGNode node, int readIndex, int readPosition);

    /** Factory preserving the original read-index-only node behavior. */
    static POGNodeFactory standard() {
        return StandardPOGNodeFactory.INSTANCE;
    }
}

final class StandardPOGNodeFactory implements POGNodeFactory {
    static final StandardPOGNodeFactory INSTANCE = new StandardPOGNodeFactory();

    private StandardPOGNodeFactory() {}

    @Override
    public POGNode create(char residue, int readIndex, int readPosition) {
        return new POGNode(residue, readIndex);
    }

    @Override
    public POGNode copy(POGNode node) {
        return new POGNode(node);
    }

    @Override
    public void addReadObservation(POGNode node, int readIndex, int readPosition) {
        node.addIDX(readIndex);
    }
}
