package jseqa.poa;

import java.util.Collections;
import java.util.List;

/** Package-private immutable snapshot of the exact dense-alignment row order. */
final class POAAlignmentTopologyView {
    final List<POGNode> nodes;
    final int[][] predecessorRows;

    POAAlignmentTopologyView(List<POGNode> nodes, int[][] predecessorRows) {
        this.nodes = Collections.unmodifiableList(nodes);
        this.predecessorRows = predecessorRows;
    }
}
