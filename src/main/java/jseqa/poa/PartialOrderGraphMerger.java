package jseqa.poa;

/** Merges two POA profiles without progressively re-adding their sequences. */
public interface PartialOrderGraphMerger {
    PartialOrderAlignment merge(PartialOrderAlignment left, PartialOrderAlignment right);
}
