package jseqa.poa;

/** One exact query-to-DAG seed anchor selected for the alignment corridor. */
public final class POASeedAnchor {
    private final int queryStart;
    private final int queryEnd;
    private final int graphStartRow;
    private final int graphEndRow;
    private final int occurrenceCount;

    POASeedAnchor(int queryStart, int queryEnd, int graphStartRow,
            int graphEndRow, int occurrenceCount) {
        this.queryStart = queryStart;
        this.queryEnd = queryEnd;
        this.graphStartRow = graphStartRow;
        this.graphEndRow = graphEndRow;
        this.occurrenceCount = occurrenceCount;
    }

    public int queryStart() { return queryStart; }
    public int queryEnd() { return queryEnd; }
    public int graphStartRow() { return graphStartRow; }
    public int graphEndRow() { return graphEndRow; }
    public int occurrenceCount() { return occurrenceCount; }
}
