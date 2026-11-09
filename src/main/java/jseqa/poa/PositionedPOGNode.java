package jseqa.poa;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalInt;

/**
 * A graph node that records the single original sequence position associated
 * with each read traversing the node.
 */
public class PositionedPOGNode extends POGNode {
    private final Map<Integer,Integer> positionsByRead;

    public PositionedPOGNode(char residue, int readIndex, int readPosition) {
        super(residue, readIndex);
        validate(readIndex, readPosition);
        this.positionsByRead = new LinkedHashMap<Integer,Integer>();
        this.positionsByRead.put(readIndex, readPosition);
    }

    PositionedPOGNode(PositionedPOGNode other) {
        super(other);
        this.positionsByRead = new LinkedHashMap<Integer,Integer>(other.positionsByRead);
    }

    /**
     * Add a read observation. Repeating the same mapping is idempotent, while
     * assigning another position from the same read to this node is invalid.
     */
    public void addReadPosition(int readIndex, int readPosition) {
        validate(readIndex, readPosition);
        Integer existing = positionsByRead.get(readIndex);
        if (existing != null) {
            if (existing.intValue() != readPosition) {
                throw new IllegalStateException("read " + readIndex
                        + " is already associated with position " + existing);
            }
            return;
        }
        positionsByRead.put(readIndex, readPosition);
        super.addIDX(readIndex);
    }

    public OptionalInt positionForRead(int readIndex) {
        Integer position = positionsByRead.get(readIndex);
        return position == null ? OptionalInt.empty() : OptionalInt.of(position);
    }

    public Map<Integer,Integer> readPositions() {
        return Collections.unmodifiableMap(positionsByRead);
    }

    @Override
    public List<Integer> getIDX() {
        return Collections.unmodifiableList(super.getIDX());
    }

    /** Position-aware nodes must be updated through {@link #addReadPosition}. */
    @Override
    public void addIDX(int readIndex) {
        throw new UnsupportedOperationException("a read position is required");
    }

    private static void validate(int readIndex, int readPosition) {
        if (readIndex < 0) {
            throw new IllegalArgumentException("read index must not be negative");
        }
        if (readPosition < 0) {
            throw new IllegalArgumentException("read position must not be negative");
        }
    }
}
