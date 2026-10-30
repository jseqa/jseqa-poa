package jseqa.poa;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/** Result of a similarity-guided agglomerative POA build. */
public final class ProgressivePOAResult {
    private final PartialOrderAlignment graph;
    private final List<String> alignment;
    private final GuideTree guideTree;
    private final List<Integer> alignmentInputIndices;

    ProgressivePOAResult(PartialOrderAlignment graph, List<String> alignment,
            GuideTree guideTree, List<Integer> alignmentInputIndices) {
        this.graph = graph;
        this.alignment = Collections.unmodifiableList(new ArrayList<String>(alignment));
        this.guideTree = guideTree;
        this.alignmentInputIndices = Collections.unmodifiableList(
                new ArrayList<Integer>(alignmentInputIndices));
    }

    public PartialOrderAlignment graph() { return graph; }
    public List<String> alignment() { return alignment; }
    public GuideTree guideTree() { return guideTree; }
    public List<Integer> alignmentInputIndices() { return alignmentInputIndices; }
}
