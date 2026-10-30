package jseqa.poa;

/**
 * Controls how a seeded aligner treats changes to its source POA graph.
 */
public enum SeededGraphMode {
    /** Check the graph before every alignment and rebuild cached data after changes. */
    MUTABLE,

    /**
     * Snapshot graph topology once when the aligner is constructed. The caller must
     * not modify the graph for the lifetime of the aligner.
     */
    IMMUTABLE
}
