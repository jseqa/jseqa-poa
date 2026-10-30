package jseqa.poa;

/** Dense fallback policy for seed-derived banded alignment. */
public enum SeededBandFallback {
    NEVER,
    DENSE_IF_NO_CHAIN,
    DENSE_IF_UNREACHABLE,
    DENSE_IF_BOUNDARY_TOUCHED
}
