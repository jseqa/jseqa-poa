package jseqa.poa;

/** Diagnostics for one seed-derived banded alignment. */
public final class SeededBandStats {
    private final int querySeeds;
    private final int repetitiveSeedsFiltered;
    private final int graphHits;
    private final int chainedAnchors;
    private final int anchoredQueryBases;
    private final int corridorGraphRows;
    private final long evaluatedCells;
    private final long denseCellCount;
    private final boolean noUsableChain;
    private final boolean endpointUnreachable;
    private final boolean corridorBoundaryTouched;
    private final boolean denseFallbackUsed;

    SeededBandStats(int querySeeds, int repetitiveSeedsFiltered, int graphHits,
            int chainedAnchors, int anchoredQueryBases, int corridorGraphRows,
            long evaluatedCells, long denseCellCount, boolean noUsableChain,
            boolean endpointUnreachable, boolean corridorBoundaryTouched,
            boolean denseFallbackUsed) {
        this.querySeeds = querySeeds;
        this.repetitiveSeedsFiltered = repetitiveSeedsFiltered;
        this.graphHits = graphHits;
        this.chainedAnchors = chainedAnchors;
        this.anchoredQueryBases = anchoredQueryBases;
        this.corridorGraphRows = corridorGraphRows;
        this.evaluatedCells = evaluatedCells;
        this.denseCellCount = denseCellCount;
        this.noUsableChain = noUsableChain;
        this.endpointUnreachable = endpointUnreachable;
        this.corridorBoundaryTouched = corridorBoundaryTouched;
        this.denseFallbackUsed = denseFallbackUsed;
    }

    public int querySeeds() { return querySeeds; }
    public int repetitiveSeedsFiltered() { return repetitiveSeedsFiltered; }
    public int graphHits() { return graphHits; }
    public int chainedAnchors() { return chainedAnchors; }
    public int anchoredQueryBases() { return anchoredQueryBases; }
    public int corridorGraphRows() { return corridorGraphRows; }
    public long evaluatedCells() { return evaluatedCells; }
    public long denseCellCount() { return denseCellCount; }
    public boolean noUsableChain() { return noUsableChain; }
    public boolean endpointUnreachable() { return endpointUnreachable; }
    public boolean corridorBoundaryTouched() { return corridorBoundaryTouched; }
    public boolean denseFallbackUsed() { return denseFallbackUsed; }
}
