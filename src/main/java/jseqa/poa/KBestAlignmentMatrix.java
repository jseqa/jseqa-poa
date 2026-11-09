package jseqa.poa;

import java.util.Arrays;
import jseqa.poa.POGPtr.PointerType;

/** Primitive storage for the k best states retained at every DP cell. */
final class KBestAlignmentMatrix {
    static final int UNREACHABLE = Integer.MIN_VALUE / 4;
    private static final PointerType[] POINTER_TYPES = PointerType.values();

    private final int rows;
    private final int columns;
    private final int ranks;
    private final int[] scores;
    private final int[] coherence;
    private final int[] previousRows;
    private final int[] previousColumns;
    private final int[] previousRanks;
    private final byte[] operations;

    KBestAlignmentMatrix(int rows, int columns, int ranks) {
        this.rows = rows;
        this.columns = columns;
        this.ranks = ranks;
        int size = rows * columns * ranks;
        scores = new int[size];
        coherence = new int[size];
        previousRows = new int[size];
        previousColumns = new int[size];
        previousRanks = new int[size];
        operations = new byte[size];
        Arrays.fill(scores, UNREACHABLE);
        Arrays.fill(previousRows, -1);
        Arrays.fill(previousColumns, -1);
        Arrays.fill(previousRanks, -1);
        Arrays.fill(operations, (byte) PointerType.UNDEF.ordinal());
    }

    int rows() {
        return rows;
    }

    int columns() {
        return columns;
    }

    int ranks() {
        return ranks;
    }

    int score(int row, int column, int rank) {
        return scores[index(row, column, rank)];
    }

    int coherence(int row, int column, int rank) {
        return coherence[index(row, column, rank)];
    }

    int previousRow(int row, int column, int rank) {
        return previousRows[index(row, column, rank)];
    }

    int previousColumn(int row, int column, int rank) {
        return previousColumns[index(row, column, rank)];
    }

    int previousRank(int row, int column, int rank) {
        return previousRanks[index(row, column, rank)];
    }

    PointerType operation(int row, int column, int rank) {
        return POINTER_TYPES[operations[index(row, column, rank)]];
    }

    void set(int row, int column, int rank, int score, int colorCoherence,
            int previousRow, int previousColumn, int previousRank, PointerType operation) {
        int index = index(row, column, rank);
        scores[index] = score;
        coherence[index] = colorCoherence;
        previousRows[index] = previousRow;
        previousColumns[index] = previousColumn;
        previousRanks[index] = previousRank;
        operations[index] = (byte) operation.ordinal();
    }

    private int index(int row, int column, int rank) {
        return ((row * columns) + column) * ranks + rank;
    }
}
