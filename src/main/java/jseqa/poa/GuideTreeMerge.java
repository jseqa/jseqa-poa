package jseqa.poa;

/** One deterministic agglomerative guide-tree merge. */
public final class GuideTreeMerge {
    private final int leftCluster;
    private final int rightCluster;
    private final int parentCluster;
    private final double similarity;
    private final int leafCount;

    GuideTreeMerge(int leftCluster, int rightCluster, int parentCluster,
            double similarity, int leafCount) {
        this.leftCluster = leftCluster;
        this.rightCluster = rightCluster;
        this.parentCluster = parentCluster;
        this.similarity = similarity;
        this.leafCount = leafCount;
    }

    public int leftCluster() { return leftCluster; }
    public int rightCluster() { return rightCluster; }
    public int parentCluster() { return parentCluster; }
    public double similarity() { return similarity; }
    public int leafCount() { return leafCount; }

    @Override
    public String toString() {
        return leftCluster + "+" + rightCluster + "->" + parentCluster
                + "@" + similarity;
    }
}
