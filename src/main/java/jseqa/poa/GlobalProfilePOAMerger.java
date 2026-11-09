package jseqa.poa;

import java.util.ArrayList;
import java.util.List;
import jseqa.alignment.substitution.SubstitutionMatrix;

/** Deterministic affine-gap global alignment and merge of two POA column profiles. */
public final class GlobalProfilePOAMerger implements PartialOrderGraphMerger {
    private static final byte MATCH = 0;
    private static final byte GAP_RIGHT = 1;
    private static final byte GAP_LEFT = 2;
    private static final double NEGATIVE_INFINITY = -Double.MAX_VALUE / 4.0;

    private final SubstitutionMatrix scoreMatrix;

    public GlobalProfilePOAMerger(SubstitutionMatrix scoreMatrix) {
        if (scoreMatrix == null) {
            throw new IllegalArgumentException("score matrix must not be null");
        }
        this.scoreMatrix = scoreMatrix;
    }

    @Override
    public PartialOrderAlignment merge(PartialOrderAlignment left,
            PartialOrderAlignment right) {
        if (left == null || right == null) {
            throw new IllegalArgumentException("both POA graphs are required");
        }
        return new POGraphGlobal(mergeAlignments(
                left.getAlignment(), right.getAlignment()), scoreMatrix);
    }

    List<String> mergeAlignments(List<String> leftRows, List<String> rightRows) {
        int leftWidth = validateAlignment(leftRows, "left");
        int rightWidth = validateAlignment(rightRows, "right");
        int columns = rightWidth + 1;
        int cells = (leftWidth + 1) * columns;
        double[] match = filled(cells, NEGATIVE_INFINITY);
        double[] gapRight = filled(cells, NEGATIVE_INFINITY);
        double[] gapLeft = filled(cells, NEGATIVE_INFINITY);
        byte[] matchTrace = new byte[cells];
        byte[] gapRightTrace = new byte[cells];
        byte[] gapLeftTrace = new byte[cells];
        match[0] = 0.0;

        for (int i = 1; i <= leftWidth; i++) {
            int current = index(i, 0, columns);
            int previous = index(i - 1, 0, columns);
            double open = match[previous] + gapCost(leftRows, i - 1,
                    scoreMatrix.gapPenalty());
            double extend = gapRight[previous] + gapCost(leftRows, i - 1,
                    scoreMatrix.gapExtensionPenalty());
            if (open >= extend) {
                gapRight[current] = open;
                gapRightTrace[current] = MATCH;
            } else {
                gapRight[current] = extend;
                gapRightTrace[current] = GAP_RIGHT;
            }
        }
        for (int j = 1; j <= rightWidth; j++) {
            int current = index(0, j, columns);
            int previous = index(0, j - 1, columns);
            double open = match[previous] + gapCost(rightRows, j - 1,
                    scoreMatrix.gapPenalty());
            double extend = gapLeft[previous] + gapCost(rightRows, j - 1,
                    scoreMatrix.gapExtensionPenalty());
            if (open >= extend) {
                gapLeft[current] = open;
                gapLeftTrace[current] = MATCH;
            } else {
                gapLeft[current] = extend;
                gapLeftTrace[current] = GAP_LEFT;
            }
        }

        for (int i = 1; i <= leftWidth; i++) {
            for (int j = 1; j <= rightWidth; j++) {
                int current = index(i, j, columns);
                int diagonal = index(i - 1, j - 1, columns);
                byte previousState = MATCH;
                double previousScore = match[diagonal];
                if (gapRight[diagonal] > previousScore) {
                    previousScore = gapRight[diagonal];
                    previousState = GAP_RIGHT;
                }
                if (gapLeft[diagonal] > previousScore) {
                    previousScore = gapLeft[diagonal];
                    previousState = GAP_LEFT;
                }
                match[current] = previousScore
                        + columnScore(leftRows, i - 1, rightRows, j - 1);
                matchTrace[current] = previousState;

                int above = index(i - 1, j, columns);
                double openRight = match[above] + gapCost(leftRows, i - 1,
                        scoreMatrix.gapPenalty());
                double extendRight = gapRight[above] + gapCost(leftRows, i - 1,
                        scoreMatrix.gapExtensionPenalty());
                double switchRight = gapLeft[above] + gapCost(leftRows, i - 1,
                        scoreMatrix.gapPenalty());
                gapRight[current] = openRight;
                gapRightTrace[current] = MATCH;
                if (extendRight > gapRight[current]) {
                    gapRight[current] = extendRight;
                    gapRightTrace[current] = GAP_RIGHT;
                }
                if (switchRight > gapRight[current]) {
                    gapRight[current] = switchRight;
                    gapRightTrace[current] = GAP_LEFT;
                }

                int before = index(i, j - 1, columns);
                double openLeft = match[before] + gapCost(rightRows, j - 1,
                        scoreMatrix.gapPenalty());
                double extendLeft = gapLeft[before] + gapCost(rightRows, j - 1,
                        scoreMatrix.gapExtensionPenalty());
                double switchLeft = gapRight[before] + gapCost(rightRows, j - 1,
                        scoreMatrix.gapPenalty());
                gapLeft[current] = openLeft;
                gapLeftTrace[current] = MATCH;
                if (extendLeft > gapLeft[current]) {
                    gapLeft[current] = extendLeft;
                    gapLeftTrace[current] = GAP_LEFT;
                }
                if (switchLeft > gapLeft[current]) {
                    gapLeft[current] = switchLeft;
                    gapLeftTrace[current] = GAP_RIGHT;
                }
            }
        }

        int end = index(leftWidth, rightWidth, columns);
        byte state = MATCH;
        double best = match[end];
        if (gapRight[end] > best) {
            best = gapRight[end];
            state = GAP_RIGHT;
        }
        if (gapLeft[end] > best) {
            state = GAP_LEFT;
        }
        return backtrack(leftRows, rightRows, leftWidth, rightWidth, columns,
                state, matchTrace, gapRightTrace, gapLeftTrace);
    }

    private List<String> backtrack(List<String> leftRows, List<String> rightRows,
            int i, int j, int columns, byte state, byte[] matchTrace,
            byte[] gapRightTrace, byte[] gapLeftTrace) {
        List<StringBuilder> merged = new ArrayList<StringBuilder>(
                leftRows.size() + rightRows.size());
        for (int row = 0; row < leftRows.size() + rightRows.size(); row++) {
            merged.add(new StringBuilder());
        }
        while (i > 0 || j > 0) {
            int current = index(i, j, columns);
            if (state == MATCH) {
                for (int row = 0; row < leftRows.size(); row++) {
                    merged.get(row).append(leftRows.get(row).charAt(i - 1));
                }
                for (int row = 0; row < rightRows.size(); row++) {
                    merged.get(leftRows.size() + row)
                            .append(rightRows.get(row).charAt(j - 1));
                }
                state = matchTrace[current];
                i--;
                j--;
            } else if (state == GAP_RIGHT) {
                for (int row = 0; row < leftRows.size(); row++) {
                    merged.get(row).append(leftRows.get(row).charAt(i - 1));
                }
                for (int row = 0; row < rightRows.size(); row++) {
                    merged.get(leftRows.size() + row).append('-');
                }
                state = gapRightTrace[current];
                i--;
            } else {
                for (int row = 0; row < leftRows.size(); row++) {
                    merged.get(row).append('-');
                }
                for (int row = 0; row < rightRows.size(); row++) {
                    merged.get(leftRows.size() + row)
                            .append(rightRows.get(row).charAt(j - 1));
                }
                state = gapLeftTrace[current];
                j--;
            }
        }
        List<String> rows = new ArrayList<String>(merged.size());
        for (StringBuilder row : merged) {
            rows.add(row.reverse().toString());
        }
        return rows;
    }

    private double columnScore(List<String> leftRows, int leftColumn,
            List<String> rightRows, int rightColumn) {
        double score = 0.0;
        for (String left : leftRows) {
            char leftBase = left.charAt(leftColumn);
            if (leftBase == '-') {
                continue;
            }
            for (String right : rightRows) {
                char rightBase = right.charAt(rightColumn);
                if (rightBase != '-') {
                    score += scoreMatrix.score(leftBase, rightBase);
                }
            }
        }
        return score / ((double) leftRows.size() * rightRows.size());
    }

    private static double gapCost(List<String> rows, int column, int penalty) {
        int residues = 0;
        for (String row : rows) {
            if (row.charAt(column) != '-') {
                residues++;
            }
        }
        return penalty * (residues / (double) rows.size());
    }

    private static int validateAlignment(List<String> rows, String side) {
        if (rows == null || rows.isEmpty()) {
            throw new IllegalArgumentException(side + " profile must contain rows");
        }
        int width = rows.get(0).length();
        if (width == 0) {
            throw new IllegalArgumentException(side + " profile rows must not be empty");
        }
        for (String row : rows) {
            if (row == null || row.length() != width) {
                throw new IllegalArgumentException(side + " profile rows must have equal width");
            }
        }
        return width;
    }

    private static double[] filled(int length, double value) {
        double[] values = new double[length];
        java.util.Arrays.fill(values, value);
        return values;
    }

    private static int index(int row, int column, int columns) {
        return row * columns + column;
    }
}
