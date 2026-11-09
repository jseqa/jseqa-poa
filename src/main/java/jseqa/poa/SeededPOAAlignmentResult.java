package jseqa.poa;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import jseqa.util.Interval;

/** Path, score, anchors, and diagnostics from a seeded banded alignment. */
public final class SeededPOAAlignmentResult {
    private final List<POGNode> path;
    private final int score;
    private final Interval queryInterval;
    private final List<POASeedAnchor> anchors;
    private final SeededBandStats stats;

    SeededPOAAlignmentResult(List<POGNode> path, int score, Interval queryInterval,
            List<POASeedAnchor> anchors, SeededBandStats stats) {
        this.path = Collections.unmodifiableList(new ArrayList<POGNode>(path));
        this.score = score;
        this.queryInterval = queryInterval;
        this.anchors = Collections.unmodifiableList(new ArrayList<POASeedAnchor>(anchors));
        this.stats = stats;
    }

    public List<POGNode> path() { return path; }
    public int score() { return score; }
    public Interval queryInterval() { return queryInterval; }
    public List<POASeedAnchor> anchors() { return anchors; }
    public SeededBandStats stats() { return stats; }
}
