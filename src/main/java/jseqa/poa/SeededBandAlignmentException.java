package jseqa.poa;

/** Indicates that a requested no-fallback seed corridor could not be aligned. */
public class SeededBandAlignmentException extends IllegalStateException {
    private static final long serialVersionUID = 1L;

    public SeededBandAlignmentException(String message) {
        super(message);
    }
}
