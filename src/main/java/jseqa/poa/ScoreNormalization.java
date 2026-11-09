package jseqa.poa;

/** Normalization applied to pairwise alignment scores before clustering. */
public enum ScoreNormalization {
    RAW,
    MIN_LENGTH,
    MAX_LENGTH,
    MEAN_LENGTH,
    SELF_SCORE
}
