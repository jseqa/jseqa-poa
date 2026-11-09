package jseqa.poa;

import java.util.Collections;
import java.util.List;
import jseqa.util.Interval;

/** One ranked sequence-to-POA alignment path. */
public final class POAAlignmentPath {
    private final List<POGNode> path;
    private final Interval queryInterval;
    private final int score;

    POAAlignmentPath(List<POGNode> path, Interval queryInterval, int score) {
        this.path = Collections.unmodifiableList(path);
        this.queryInterval = queryInterval;
        this.score = score;
    }

    /** @return graph nodes aligned to successive query positions */
    public List<POGNode> path() {
        return path;
    }

    /** @return zero-based, half-open aligned query interval */
    public Interval queryInterval() {
        return queryInterval;
    }

    /** @return primary sequence-to-graph alignment score */
    public int score() {
        return score;
    }
}
