package jseqa.poa;

import jseqa.alignment.AlignmentMatrixScore;

/** Alignment score matrix with reusable grow-only backing storage. */
final class ReusableAlignmentMatrixScore extends AlignmentMatrixScore {
    private int[] values;
    private int capacityRows;
    private int capacityColumns;

    ReusableAlignmentMatrixScore(int rows, int columns) {
        super(1, 1);
        resize(rows, columns);
    }

    void resize(int rows, int columns) {
        if (rows < 1 || columns < 1) {
            throw new IllegalArgumentException("matrix dimensions must be positive");
        }
        if (values == null || rows > capacityRows || columns > capacityColumns) {
            int newRows = grow(capacityRows, rows);
            int newColumns = grow(capacityColumns, columns);
            values = new int[newRows * newColumns];
            capacityRows = newRows;
            capacityColumns = newColumns;
        }
        this.numRows = rows;
        this.numCols = columns;
    }

    @Override
    public int get(int row, int column) {
        return values[index(row, column)];
    }

    @Override
    public void set(int row, int column, int value) {
        values[index(row, column)] = value;
    }

    @Override
    public void clearCell(int row, int column) {
        // AlignmentMatrixScore's constructor invokes this virtual method before
        // this subclass has initialized its reusable backing array.
        if (values != null) {
            values[index(row, column)] = 0;
        }
    }

    @Override
    public void setLowestVal(int row, int column) {
        values[index(row, column)] = Integer.MIN_VALUE;
    }

    @Override
    public String cellToString(int row, int column) {
        return Integer.toString(get(row, column));
    }

    private int index(int row, int column) {
        return row * capacityColumns + column;
    }

    private static int grow(int current, int required) {
        if (current >= required) {
            return current;
        }
        int grown = current == 0 ? 1 : current;
        while (grown < required) {
            grown = grown < 1024 ? grown * 2 : grown + grown / 2;
        }
        return grown;
    }
}
