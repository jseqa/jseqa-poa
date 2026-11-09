package jseqa.poa;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/** Immutable, inspectable guide tree for an agglomerative POA build. */
public final class GuideTree {
    /** Immutable guide-tree node. */
    public static final class Node {
        private final int id;
        private final int inputIndex;
        private final Node left;
        private final Node right;
        private final List<Integer> leafIndices;
        private final double mergeSimilarity;

        static Node leaf(int inputIndex) {
            return new Node(inputIndex, inputIndex, null, null,
                    Collections.singletonList(inputIndex), Double.NaN);
        }

        static Node parent(int id, Node left, Node right, double similarity) {
            List<Integer> leaves = new ArrayList<Integer>(
                    left.leafIndices.size() + right.leafIndices.size());
            leaves.addAll(left.leafIndices);
            leaves.addAll(right.leafIndices);
            Collections.sort(leaves);
            return new Node(id, -1, left, right,
                    Collections.unmodifiableList(leaves), similarity);
        }

        private Node(int id, int inputIndex, Node left, Node right,
                List<Integer> leaves, double similarity) {
            this.id = id;
            this.inputIndex = inputIndex;
            this.left = left;
            this.right = right;
            this.leafIndices = leaves;
            this.mergeSimilarity = similarity;
        }

        public int id() { return id; }
        public boolean isLeaf() { return inputIndex >= 0; }
        public int inputIndex() { return inputIndex; }
        public Node left() { return left; }
        public Node right() { return right; }
        public List<Integer> leafIndices() { return leafIndices; }
        public int leafCount() { return leafIndices.size(); }
        public double mergeSimilarity() { return mergeSimilarity; }
    }

    private final Node root;
    private final List<GuideTreeMerge> merges;

    GuideTree(Node root, List<GuideTreeMerge> merges) {
        this.root = root;
        this.merges = Collections.unmodifiableList(
                new ArrayList<GuideTreeMerge>(merges));
    }

    public Node root() { return root; }
    public List<GuideTreeMerge> merges() { return merges; }

    public String toNewick() {
        return newick(root) + ";";
    }

    private static String newick(Node node) {
        if (node.isLeaf()) {
            return Integer.toString(node.inputIndex());
        }
        return "(" + newick(node.left()) + "," + newick(node.right()) + ")"
                + node.id();
    }
}
