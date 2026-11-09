package jseqa.poa;

/** Immutable configuration for a seed-derived DAG alignment corridor. */
public final class SeededBandOptions {
    private final int kmerLength;
    private final SeedSelection seedSelection;
    private final int minimizerWindow;
    private final int queryPadding;
    private final int graphPadding;
    private final int maxSeedOccurrences;
    private final int maxSeedHits;
    private final int minimumChainedSeeds;
    private final SeededBandFallback fallback;

    private SeededBandOptions(Builder builder) {
        this.kmerLength = builder.kmerLength;
        this.seedSelection = builder.seedSelection;
        this.minimizerWindow = builder.minimizerWindow;
        this.queryPadding = builder.queryPadding;
        this.graphPadding = builder.graphPadding;
        this.maxSeedOccurrences = builder.maxSeedOccurrences;
        this.maxSeedHits = builder.maxSeedHits;
        this.minimumChainedSeeds = builder.minimumChainedSeeds;
        this.fallback = builder.fallback;
    }

    public static Builder builder() { return new Builder(); }
    public int kmerLength() { return kmerLength; }
    public SeedSelection seedSelection() { return seedSelection; }
    public int minimizerWindow() { return minimizerWindow; }
    public int queryPadding() { return queryPadding; }
    public int graphPadding() { return graphPadding; }
    public int maxSeedOccurrences() { return maxSeedOccurrences; }
    public int maxSeedHits() { return maxSeedHits; }
    public int minimumChainedSeeds() { return minimumChainedSeeds; }
    public SeededBandFallback fallback() { return fallback; }

    public static final class Builder {
        private int kmerLength = 11;
        private SeedSelection seedSelection = SeedSelection.MINIMIZERS;
        private int minimizerWindow = 8;
        private int queryPadding = 24;
        private int graphPadding = 8;
        private int maxSeedOccurrences = 64;
        private int maxSeedHits = 2000;
        private int minimumChainedSeeds = 2;
        private SeededBandFallback fallback = SeededBandFallback.DENSE_IF_UNREACHABLE;

        public Builder kmerLength(int value) { this.kmerLength = value; return this; }
        public Builder seedSelection(SeedSelection value) { this.seedSelection = value; return this; }
        public Builder minimizerWindow(int value) { this.minimizerWindow = value; return this; }
        public Builder queryPadding(int value) { this.queryPadding = value; return this; }
        public Builder graphPadding(int value) { this.graphPadding = value; return this; }
        public Builder maxSeedOccurrences(int value) { this.maxSeedOccurrences = value; return this; }
        public Builder maxSeedHits(int value) { this.maxSeedHits = value; return this; }
        public Builder minimumChainedSeeds(int value) { this.minimumChainedSeeds = value; return this; }
        public Builder fallback(SeededBandFallback value) { this.fallback = value; return this; }

        public SeededBandOptions build() {
            if (kmerLength < 2 || minimizerWindow < 1 || queryPadding < 0
                    || graphPadding < 0 || maxSeedOccurrences < 1
                    || maxSeedHits < 1 || minimumChainedSeeds < 1
                    || seedSelection == null || fallback == null) {
                throw new IllegalArgumentException("invalid seeded-band configuration");
            }
            return new SeededBandOptions(this);
        }
    }
}
