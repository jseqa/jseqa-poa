package jseqa.poa;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import jseqa.alignment.substitution.SubstitutionMatrix;

/** Similarity-guided, profile-to-profile agglomerative POA builder. */
public final class AgglomerativePOABuilder {
    private static final int NEGATIVE_INFINITY = Integer.MIN_VALUE / 4;

    private final SubstitutionMatrix scoreMatrix;
    private final Linkage linkage;
    private final ScoreNormalization normalization;
    private PairwiseAlignmentMode pairwiseMode = PairwiseAlignmentMode.GLOBAL;
    private int threads = 1;

    public AgglomerativePOABuilder(SubstitutionMatrix scoreMatrix) {
        this(scoreMatrix, Linkage.UPGMA, ScoreNormalization.SELF_SCORE);
    }

    public AgglomerativePOABuilder(SubstitutionMatrix scoreMatrix,
            Linkage linkage, ScoreNormalization normalization) {
        if (scoreMatrix == null || linkage == null || normalization == null) {
            throw new IllegalArgumentException(
                    "score matrix, linkage, and normalization are required");
        }
        this.scoreMatrix = scoreMatrix;
        this.linkage = linkage;
        this.normalization = normalization;
    }

    public AgglomerativePOABuilder pairwiseMode(PairwiseAlignmentMode mode) {
        if (mode == null) {
            throw new IllegalArgumentException("pairwise alignment mode is required");
        }
        this.pairwiseMode = mode;
        return this;
    }

    public AgglomerativePOABuilder threads(int threadCount) {
        if (threadCount < 1) {
            throw new IllegalArgumentException("thread count must be positive");
        }
        this.threads = threadCount;
        return this;
    }

    public ProgressivePOAResult build(List<String> sequences) {
        List<String> inputs = validateSequences(sequences);
        PairwiseScores pairwise = computePairwiseScores(inputs);
        GuideTree tree = buildGuideTree(inputs.size(), pairwise);
        Profile profile = buildProfile(tree.root(), inputs,
                new GlobalProfilePOAMerger(scoreMatrix));
        List<String> alignment = profile.graph.getAlignment();
        validateResult(inputs, alignment, profile.inputIndices);
        return new ProgressivePOAResult(profile.graph, alignment, tree,
                profile.inputIndices);
    }

    /** Compute the symmetric normalized similarity matrix used for clustering. */
    public double[][] pairwiseSimilarities(List<String> sequences) {
        List<String> inputs = validateSequences(sequences);
        PairwiseScores scores = computePairwiseScores(inputs);
        double[][] matrix = new double[inputs.size()][inputs.size()];
        for (int i = 0; i < inputs.size(); i++) {
            matrix[i][i] = normalize(scores.selfScores[i], scores.selfScores[i],
                    scores.selfScores[i], inputs.get(i).length(), inputs.get(i).length());
            for (int j = i + 1; j < inputs.size(); j++) {
                matrix[i][j] = scores.get(i, j);
                matrix[j][i] = matrix[i][j];
            }
        }
        return matrix;
    }

    public static ProgressivePOAResult agglomerative(List<String> sequences,
            SubstitutionMatrix scoreMatrix) {
        return new AgglomerativePOABuilder(scoreMatrix).build(sequences);
    }

    private PairwiseScores computePairwiseScores(final List<String> sequences) {
        final int size = sequences.size();
        final int[] selfScores = new int[size];
        final double[] similarities = new double[size * (size - 1) / 2];
        ExecutorService executor = threads == 1 ? null
                : Executors.newFixedThreadPool(threads);
        try {
            if (executor == null) {
                for (int i = 0; i < size; i++) {
                    selfScores[i] = alignmentScore(sequences.get(i), sequences.get(i),
                            PairwiseAlignmentMode.GLOBAL);
                }
            } else {
                List<Callable<Integer>> selfTasks = new ArrayList<Callable<Integer>>(size);
                for (int i = 0; i < size; i++) {
                    final String sequence = sequences.get(i);
                    selfTasks.add(new Callable<Integer>() {
                        @Override public Integer call() {
                            return alignmentScore(sequence, sequence,
                                    PairwiseAlignmentMode.GLOBAL);
                        }
                    });
                }
                List<Future<Integer>> futures = executor.invokeAll(selfTasks);
                for (int i = 0; i < futures.size(); i++) {
                    selfScores[i] = futures.get(i).get();
                }
            }

            List<Callable<Void>> tasks = new ArrayList<Callable<Void>>(
                    size * (size - 1) / 2);
            for (int i = 0; i < size; i++) {
                for (int j = i + 1; j < size; j++) {
                    final int left = i;
                    final int right = j;
                    Callable<Void> task = new Callable<Void>() {
                        @Override public Void call() {
                            int score;
                            if (pairwiseMode == PairwiseAlignmentMode.FITTING) {
                                score = Math.max(
                                        alignmentScore(sequences.get(left), sequences.get(right),
                                                PairwiseAlignmentMode.FITTING),
                                        alignmentScore(sequences.get(right), sequences.get(left),
                                                PairwiseAlignmentMode.FITTING));
                            } else {
                                score = alignmentScore(sequences.get(left), sequences.get(right),
                                        pairwiseMode);
                            }
                            similarities[triangleIndex(left, right, size)] = normalize(score,
                                    selfScores[left], selfScores[right],
                                    sequences.get(left).length(), sequences.get(right).length());
                            return null;
                        }
                    };
                    if (executor == null) {
                        try {
                            task.call();
                        } catch (Exception impossible) {
                            throw new IllegalStateException(impossible);
                        }
                    } else {
                        tasks.add(task);
                    }
                }
            }
            if (executor != null) {
                for (Future<Void> future : executor.invokeAll(tasks)) {
                    future.get();
                }
            }
            return new PairwiseScores(size, similarities, selfScores);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IllegalStateException("pairwise scoring was interrupted", e);
        } catch (Exception e) {
            throw new IllegalStateException("pairwise scoring failed", e);
        } finally {
            if (executor != null) {
                executor.shutdown();
            }
        }
    }

    private int alignmentScore(String left, String right, PairwiseAlignmentMode mode) {
        int rightLength = right.length();
        int[] previousMatch = filled(rightLength + 1, NEGATIVE_INFINITY);
        int[] previousGapRight = filled(rightLength + 1, NEGATIVE_INFINITY);
        int[] previousGapLeft = filled(rightLength + 1, NEGATIVE_INFINITY);
        previousMatch[0] = 0;
        boolean freeLeftBoundary = mode == PairwiseAlignmentMode.OVERLAP
                || mode == PairwiseAlignmentMode.FITTING
                || mode == PairwiseAlignmentMode.LOCAL;
        boolean freeTopBoundary = mode == PairwiseAlignmentMode.OVERLAP
                || mode == PairwiseAlignmentMode.LOCAL;
        for (int j = 1; j <= rightLength; j++) {
            if (freeTopBoundary) {
                previousMatch[j] = 0;
            } else {
                previousGapLeft[j] = scoreMatrix.gapPenalty()
                        + (j - 1) * scoreMatrix.gapExtensionPenalty();
            }
        }
        int bestLocal = 0;
        int bestLastColumn = max(previousMatch[rightLength],
                previousGapRight[rightLength], previousGapLeft[rightLength]);
        int[] currentMatch = new int[rightLength + 1];
        int[] currentGapRight = new int[rightLength + 1];
        int[] currentGapLeft = new int[rightLength + 1];
        for (int i = 1; i <= left.length(); i++) {
            Arrays.fill(currentMatch, NEGATIVE_INFINITY);
            Arrays.fill(currentGapRight, NEGATIVE_INFINITY);
            Arrays.fill(currentGapLeft, NEGATIVE_INFINITY);
            if (freeLeftBoundary) {
                currentMatch[0] = 0;
            } else {
                currentGapRight[0] = scoreMatrix.gapPenalty()
                        + (i - 1) * scoreMatrix.gapExtensionPenalty();
            }
            for (int j = 1; j <= rightLength; j++) {
                currentMatch[j] = max(previousMatch[j - 1],
                        previousGapRight[j - 1], previousGapLeft[j - 1])
                        + scoreMatrix.score(left.charAt(i - 1), right.charAt(j - 1));
                currentGapRight[j] = max(
                        previousMatch[j] + scoreMatrix.gapPenalty(),
                        previousGapRight[j] + scoreMatrix.gapExtensionPenalty(),
                        previousGapLeft[j] + scoreMatrix.gapPenalty());
                currentGapLeft[j] = max(
                        currentMatch[j - 1] + scoreMatrix.gapPenalty(),
                        currentGapLeft[j - 1] + scoreMatrix.gapExtensionPenalty(),
                        currentGapRight[j - 1] + scoreMatrix.gapPenalty());
                if (mode == PairwiseAlignmentMode.LOCAL) {
                    currentMatch[j] = Math.max(0, currentMatch[j]);
                    currentGapRight[j] = Math.max(0, currentGapRight[j]);
                    currentGapLeft[j] = Math.max(0, currentGapLeft[j]);
                    bestLocal = Math.max(bestLocal, max(currentMatch[j],
                            currentGapRight[j], currentGapLeft[j]));
                }
            }
            bestLastColumn = Math.max(bestLastColumn, max(currentMatch[rightLength],
                    currentGapRight[rightLength], currentGapLeft[rightLength]));
            int[] swap = previousMatch; previousMatch = currentMatch; currentMatch = swap;
            swap = previousGapRight; previousGapRight = currentGapRight; currentGapRight = swap;
            swap = previousGapLeft; previousGapLeft = currentGapLeft; currentGapLeft = swap;
        }
        if (mode == PairwiseAlignmentMode.LOCAL) {
            return bestLocal;
        }
        if (mode == PairwiseAlignmentMode.FITTING) {
            return bestLastColumn;
        }
        if (mode == PairwiseAlignmentMode.OVERLAP) {
            int bestLastRow = NEGATIVE_INFINITY;
            for (int j = 0; j <= rightLength; j++) {
                bestLastRow = Math.max(bestLastRow, max(previousMatch[j],
                        previousGapRight[j], previousGapLeft[j]));
            }
            return Math.max(bestLastColumn, bestLastRow);
        }
        return max(previousMatch[rightLength], previousGapRight[rightLength],
                previousGapLeft[rightLength]);
    }

    private double normalize(int score, int leftSelf, int rightSelf,
            int leftLength, int rightLength) {
        switch (normalization) {
        case RAW:
            return score;
        case MIN_LENGTH:
            return score / (double) Math.min(leftLength, rightLength);
        case MAX_LENGTH:
            return score / (double) Math.max(leftLength, rightLength);
        case MEAN_LENGTH:
            return score / ((leftLength + rightLength) / 2.0);
        case SELF_SCORE:
        default:
            int denominator = leftSelf + rightSelf;
            return denominator == 0 ? (score == 0 ? 1.0 : 0.0)
                    : (2.0 * score) / denominator;
        }
    }

    private GuideTree buildGuideTree(int sequenceCount, PairwiseScores similarities) {
        List<Cluster> active = new ArrayList<Cluster>(sequenceCount);
        for (int i = 0; i < sequenceCount; i++) {
            active.add(new Cluster(i, GuideTree.Node.leaf(i),
                    Collections.singletonList(i)));
        }
        List<GuideTreeMerge> merges = new ArrayList<GuideTreeMerge>();
        int nextId = sequenceCount;
        while (active.size() > 1) {
            Cluster bestLeft = null;
            Cluster bestRight = null;
            double bestSimilarity = -Double.MAX_VALUE;
            for (int i = 0; i < active.size(); i++) {
                for (int j = i + 1; j < active.size(); j++) {
                    Cluster left = canonicalFirst(active.get(i), active.get(j));
                    Cluster right = left == active.get(i) ? active.get(j) : active.get(i);
                    double similarity = clusterSimilarity(left, right, similarities);
                    if (bestLeft == null || similarity > bestSimilarity
                            || (Double.compare(similarity, bestSimilarity) == 0
                            && pairComesFirst(left, right, bestLeft, bestRight))) {
                        bestLeft = left;
                        bestRight = right;
                        bestSimilarity = similarity;
                    }
                }
            }
            GuideTree.Node parentNode = GuideTree.Node.parent(nextId,
                    bestLeft.node, bestRight.node, bestSimilarity);
            List<Integer> leaves = new ArrayList<Integer>(bestLeft.leaves.size()
                    + bestRight.leaves.size());
            leaves.addAll(bestLeft.leaves);
            leaves.addAll(bestRight.leaves);
            Collections.sort(leaves);
            Cluster parent = new Cluster(nextId, parentNode,
                    Collections.unmodifiableList(leaves));
            merges.add(new GuideTreeMerge(bestLeft.id, bestRight.id, nextId,
                    bestSimilarity, leaves.size()));
            active.remove(bestLeft);
            active.remove(bestRight);
            active.add(parent);
            nextId++;
        }
        return new GuideTree(active.get(0).node, merges);
    }

    private double clusterSimilarity(Cluster left, Cluster right,
            PairwiseScores similarities) {
        double result = linkage == Linkage.COMPLETE ? Double.MAX_VALUE
                : (linkage == Linkage.SINGLE ? -Double.MAX_VALUE : 0.0);
        int pairs = 0;
        for (int leftLeaf : left.leaves) {
            for (int rightLeaf : right.leaves) {
                double value = similarities.get(leftLeaf, rightLeaf);
                if (linkage == Linkage.SINGLE) {
                    result = Math.max(result, value);
                } else if (linkage == Linkage.COMPLETE) {
                    result = Math.min(result, value);
                } else {
                    result += value;
                    pairs++;
                }
            }
        }
        return linkage == Linkage.UPGMA || linkage == Linkage.AVERAGE
                ? result / pairs : result;
    }

    private Profile buildProfile(GuideTree.Node node, List<String> sequences,
            GlobalProfilePOAMerger merger) {
        if (node.isLeaf()) {
            return new Profile(new POGraphGlobal(
                    sequences.get(node.inputIndex()), scoreMatrix),
                    Collections.singletonList(node.inputIndex()));
        }
        Profile left = buildProfile(node.left(), sequences, merger);
        Profile right = buildProfile(node.right(), sequences, merger);
        PartialOrderAlignment merged = merger.merge(left.graph, right.graph);
        List<Integer> indices = new ArrayList<Integer>(
                left.inputIndices.size() + right.inputIndices.size());
        indices.addAll(left.inputIndices);
        indices.addAll(right.inputIndices);
        validateResult(sequences, merged.getAlignment(), indices);
        return new Profile(merged, Collections.unmodifiableList(indices));
    }

    private static void validateResult(List<String> inputs, List<String> alignment,
            List<Integer> inputIndices) {
        if (alignment.size() != inputIndices.size()) {
            throw new IllegalStateException("agglomerative merge lost or duplicated rows");
        }
        boolean[] seen = new boolean[inputs.size()];
        int width = alignment.get(0).length();
        for (int row = 0; row < alignment.size(); row++) {
            int input = inputIndices.get(row);
            if (input < 0 || input >= inputs.size() || seen[input]) {
                throw new IllegalStateException("invalid or duplicate input-row mapping");
            }
            seen[input] = true;
            if (alignment.get(row).length() != width
                    || !alignment.get(row).replace("-", "").equals(inputs.get(input))) {
                throw new IllegalStateException("merged row does not preserve input " + input);
            }
        }
    }

    private static List<String> validateSequences(List<String> sequences) {
        if (sequences == null || sequences.isEmpty()) {
            throw new IllegalArgumentException("at least one sequence is required");
        }
        List<String> copy = new ArrayList<String>(sequences.size());
        for (String sequence : sequences) {
            if (sequence == null || sequence.isEmpty()) {
                throw new IllegalArgumentException("sequences must not be null or empty");
            }
            for (int i = 0; i < sequence.length(); i++) {
                if (!Character.isLetter(sequence.charAt(i))) {
                    throw new IllegalArgumentException(
                            "invalid sequence symbol: " + sequence.charAt(i));
                }
            }
            copy.add(sequence);
        }
        return Collections.unmodifiableList(copy);
    }

    private static Cluster canonicalFirst(Cluster left, Cluster right) {
        int comparison = Integer.compare(left.minLeaf(), right.minLeaf());
        if (comparison != 0) {
            return comparison < 0 ? left : right;
        }
        return left.id < right.id ? left : right;
    }

    private static boolean pairComesFirst(Cluster left, Cluster right,
            Cluster currentLeft, Cluster currentRight) {
        int comparison = Integer.compare(left.minLeaf(), currentLeft.minLeaf());
        if (comparison != 0) return comparison < 0;
        comparison = Integer.compare(right.minLeaf(), currentRight.minLeaf());
        if (comparison != 0) return comparison < 0;
        comparison = Integer.compare(left.id, currentLeft.id);
        if (comparison != 0) return comparison < 0;
        return right.id < currentRight.id;
    }

    private static int triangleIndex(int left, int right, int size) {
        if (left > right) {
            int swap = left; left = right; right = swap;
        }
        return left * size - left * (left + 1) / 2 + right - left - 1;
    }

    private static int[] filled(int length, int value) {
        int[] values = new int[length];
        Arrays.fill(values, value);
        return values;
    }

    private static int max(int first, int second, int third) {
        return Math.max(first, Math.max(second, third));
    }

    private static final class PairwiseScores {
        final int size;
        final double[] values;
        final int[] selfScores;

        PairwiseScores(int size, double[] values, int[] selfScores) {
            this.size = size;
            this.values = values;
            this.selfScores = selfScores;
        }

        double get(int left, int right) {
            if (left == right) return 1.0;
            return values[triangleIndex(left, right, size)];
        }
    }

    private static final class Cluster {
        final int id;
        final GuideTree.Node node;
        final List<Integer> leaves;

        Cluster(int id, GuideTree.Node node, List<Integer> leaves) {
            this.id = id;
            this.node = node;
            this.leaves = leaves;
        }

        int minLeaf() { return leaves.get(0); }
    }

    private static final class Profile {
        final PartialOrderAlignment graph;
        final List<Integer> inputIndices;

        Profile(PartialOrderAlignment graph, List<Integer> inputIndices) {
            this.graph = graph;
            this.inputIndices = inputIndices;
        }
    }
}
