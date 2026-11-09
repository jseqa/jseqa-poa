package jseqa.poa;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import jseqa.alignment.substitution.SubstitutionMatrix;
import jseqa.poa.POGPtr.PointerType;
import jseqa.util.Interval;
import org.jgrapht.graph.DefaultDirectedGraph;

/**
 * Independent seed-derived banded global or local aligner for a sequence and an
 * existing POA DAG. It does not modify the POA graph or its dense alignment
 * workspaces.
 */
public final class SeededBandedPOAAligner {
    private static final int UNREACHABLE = Integer.MIN_VALUE / 4;
    private static final PointerType[] POINTER_TYPES = PointerType.values();

    private final PartialOrderAlignment poa;
    private final SubstitutionMatrix scores;
    private final boolean localMode;
    private final SeededGraphMode graphMode;
    private Topology topology;
    private long topologyFingerprint = Long.MIN_VALUE;
    private final Map<IndexKey,SeedIndex> indexCache =
            new HashMap<IndexKey,SeedIndex>();
    private final ThreadLocal<Workspace> workspaces =
            new ThreadLocal<Workspace>();

    public SeededBandedPOAAligner(PartialOrderAlignment poa,
            SubstitutionMatrix scores) {
        this(poa, scores, SeededGraphMode.MUTABLE);
    }

    /**
     * Creates an aligner with explicit graph cache semantics.
     *
     * <p>{@link SeededGraphMode#IMMUTABLE} avoids the per-query graph
     * fingerprint traversal. In that mode the graph must not be modified after
     * this constructor returns.</p>
     */
    public SeededBandedPOAAligner(PartialOrderAlignment poa,
            SubstitutionMatrix scores, SeededGraphMode graphMode) {
        if (poa == null || scores == null) {
            throw new IllegalArgumentException("POA graph and score matrix are required");
        }
        if (graphMode == null) {
            throw new IllegalArgumentException("seeded graph mode is required");
        }
        if (!(poa instanceof POGraphGlobal) && !(poa instanceof POGraphLocal)) {
            throw new UnsupportedOperationException(
                    "seed-derived banding currently supports global and local POA graphs");
        }
        this.poa = poa;
        this.scores = scores;
        this.localMode = poa instanceof POGraphLocal;
        this.graphMode = graphMode;
        if (graphMode == SeededGraphMode.IMMUTABLE) {
            topology = buildTopology(poa.getGraph());
        }
    }

    public SeededPOAAlignmentResult align(String query, SeededBandOptions options) {
        if (query == null || query.isEmpty() || options == null) {
            throw new IllegalArgumentException("query and seeded-band options are required");
        }
        refreshTopology();
        if (topology.nodes.isEmpty()) {
            throw new IllegalArgumentException("POA graph must contain at least one node");
        }
        SeedIndex index = seedIndex(options);
        SeedSearch search = findAnchors(query, options, index);
        long denseCells = (long) topology.nodes.size() * query.length();
        if (search.chain.size() < options.minimumChainedSeeds()) {
            if (options.fallback() == SeededBandFallback.NEVER) {
                throw new SeededBandAlignmentException("no usable seed chain");
            }
            AlignmentRun dense = runExistingDense(query, 0);
            return result(dense, search, denseCells, true, true, false,
                    topology.nodes.size());
        }

        RowIntervals corridor = buildCorridor(query.length(), search.chain, options);
        AlignmentRun banded = run(query, corridor);
        boolean fallback = false;
        if (!banded.reachable && (options.fallback() == SeededBandFallback.DENSE_IF_UNREACHABLE
                || options.fallback() == SeededBandFallback.DENSE_IF_BOUNDARY_TOUCHED)) {
            fallback = true;
        } else if (banded.boundaryTouched
                && options.fallback() == SeededBandFallback.DENSE_IF_BOUNDARY_TOUCHED) {
            fallback = true;
        }
        if (!banded.reachable && !fallback) {
            throw new SeededBandAlignmentException("seed-derived corridor has no global endpoint");
        }
        if (fallback) {
            AlignmentRun dense = runExistingDense(query, banded.evaluatedCells);
            return result(dense, search, denseCells, true, false,
                    !banded.reachable, corridor.activeRows());
        }
        return result(banded, search, denseCells, false, false, false,
                corridor.activeRows());
    }

    public int cachedSeedIndexCount() {
        return indexCache.size();
    }

    private SeededPOAAlignmentResult result(AlignmentRun run, SeedSearch search,
            long denseCells, boolean fallback, boolean noChain,
            boolean unreachable, int corridorRows) {
        int anchoredBases = 0;
        if (!search.chain.isEmpty()) {
            anchoredBases = search.chain.get(search.chain.size() - 1).queryEnd()
                    - search.chain.get(0).queryStart();
        }
        SeededBandStats stats = new SeededBandStats(search.querySeeds,
                search.repetitiveFiltered, search.hitCount, search.chain.size(),
                anchoredBases, corridorRows, run.evaluatedCells, denseCells,
                noChain, unreachable, run.boundaryTouched, fallback);
        return new SeededPOAAlignmentResult(run.path, run.score,
                run.queryInterval, search.chain, stats);
    }

    private void refreshTopology() {
        if (graphMode == SeededGraphMode.IMMUTABLE) return;
        DefaultDirectedGraph<POGNode,POGEdge> graph = poa.getGraph();
        long fingerprint = fingerprint(graph);
        if (topology == null || fingerprint != topologyFingerprint) {
            topology = buildTopology(graph);
            topologyFingerprint = fingerprint;
            indexCache.clear();
        }
    }

    private Topology buildTopology(DefaultDirectedGraph<POGNode,POGEdge> graph) {
        return Topology.build(graph,
                ((AbstractPartialOrderGraph) poa).alignmentTopologyView());
    }

    private SeedIndex seedIndex(SeededBandOptions options) {
        IndexKey key = new IndexKey(options.kmerLength(), options.seedSelection(),
                options.minimizerWindow());
        SeedIndex index = indexCache.get(key);
        if (index == null) {
            index = SeedIndex.build(topology, options);
            indexCache.put(key, index);
        }
        return index;
    }

    private SeedSearch findAnchors(String query, SeededBandOptions options,
            SeedIndex index) {
        List<SeedToken> querySeeds = selectSeeds(query, options.kmerLength(),
                options.seedSelection(), options.minimizerWindow());
        List<Hit> hits = new ArrayList<Hit>();
        int repetitive = 0;
        for (SeedToken token : querySeeds) {
            List<Occurrence> occurrences = index.byKmer.get(token.kmer);
            if (occurrences == null) continue;
            if (occurrences.size() > options.maxSeedOccurrences()) {
                repetitive++;
                continue;
            }
            for (Occurrence occurrence : occurrences) {
                hits.add(new Hit(token.position, occurrence, occurrences.size(),
                        options.kmerLength()));
            }
        }
        Collections.sort(hits, Hit.RARITY_ORDER);
        if (hits.size() > options.maxSeedHits()) {
            hits = new ArrayList<Hit>(hits.subList(0, options.maxSeedHits()));
        }
        Collections.sort(hits, Hit.CHAIN_ORDER);
        List<POASeedAnchor> chain = chain(hits, options.kmerLength());
        return new SeedSearch(querySeeds.size(), repetitive, hits.size(), chain);
    }

    private List<POASeedAnchor> chain(List<Hit> hits, int k) {
        if (hits.isEmpty()) return Collections.emptyList();
        int[] best = new int[hits.size()];
        int[] previous = new int[hits.size()];
        int[] count = new int[hits.size()];
        Arrays.fill(previous, -1);
        int bestIndex = 0;
        for (int i = 0; i < hits.size(); i++) {
            Hit current = hits.get(i);
            best[i] = k * 1000 / current.occurrenceCount;
            count[i] = 1;
            for (int j = 0; j < i; j++) {
                Hit prior = hits.get(j);
                if (prior.queryStart >= current.queryStart
                        || !topology.reaches(prior.occurrence.startRow,
                                current.occurrence.startRow)) {
                    continue;
                }
                int queryDelta = current.queryStart - prior.queryStart;
                int graphDelta = topology.forwardDepth[current.occurrence.startRow - 1]
                        - topology.forwardDepth[prior.occurrence.startRow - 1];
                int candidate = best[j] + k * 1000 / current.occurrenceCount
                        - Math.abs(queryDelta - graphDelta) * 10;
                if (candidate > best[i]
                        || (candidate == best[i] && count[j] + 1 > count[i])) {
                    best[i] = candidate;
                    previous[i] = j;
                    count[i] = count[j] + 1;
                }
            }
            if (best[i] > best[bestIndex]
                    || (best[i] == best[bestIndex] && count[i] > count[bestIndex])) {
                bestIndex = i;
            }
        }
        List<POASeedAnchor> anchors = new ArrayList<POASeedAnchor>();
        for (int index = bestIndex; index >= 0; index = previous[index]) {
            Hit hit = hits.get(index);
            anchors.add(new POASeedAnchor(hit.queryStart, hit.queryStart + k,
                    hit.occurrence.startRow, hit.occurrence.endRow,
                    hit.occurrenceCount));
            if (previous[index] < 0) break;
        }
        Collections.reverse(anchors);
        return anchors;
    }

    private RowIntervals buildCorridor(int queryLength, List<POASeedAnchor> anchors,
            SeededBandOptions options) {
        int rows = topology.nodes.size();
        boolean[] included = new boolean[rows + 1];
        POASeedAnchor first = anchors.get(0);
        POASeedAnchor last = anchors.get(anchors.size() - 1);
        for (int row = 1; row <= rows; row++) {
            if (!localMode && (topology.reaches(row, first.graphStartRow())
                    || topology.reaches(last.graphStartRow(), row))) {
                included[row] = true;
            }
            for (int anchor = 0; anchor + 1 < anchors.size(); anchor++) {
                if (topology.reaches(anchors.get(anchor).graphStartRow(), row)
                        && topology.reaches(row,
                                anchors.get(anchor + 1).graphStartRow())) {
                    included[row] = true;
                    break;
                }
            }
        }
        for (POASeedAnchor anchor : anchors) {
            for (int row = anchor.graphStartRow(); row <= anchor.graphEndRow(); row++) {
                included[row] = true;
            }
        }
        for (int expansion = 0; expansion < options.graphPadding(); expansion++) {
            boolean[] expanded = included.clone();
            for (int row = 1; row <= rows; row++) {
                if (!included[row]) continue;
                for (int predecessor : topology.predecessors[row - 1]) expanded[predecessor] = true;
                for (int successor : topology.successors[row - 1]) expanded[successor] = true;
            }
            included = expanded;
        }

        int[] starts = new int[rows + 1];
        int[] ends = new int[rows + 1];
        Arrays.fill(starts, 1);
        for (int row = 1; row <= rows; row++) {
            if (!included[row]) continue;
            int center = projectedQueryPosition(row, queryLength, anchors);
            starts[row] = Math.max(1, center - options.queryPadding());
            ends[row] = Math.min(queryLength, center + options.queryPadding());
            if (!localMode && topology.predecessors[row - 1].length == 0) starts[row] = 1;
            if (!localMode && topology.successors[row - 1].length == 0) ends[row] = queryLength;
        }
        return new RowIntervals(starts, ends);
    }

    private int projectedQueryPosition(int row, int queryLength,
            List<POASeedAnchor> anchors) {
        int depth = topology.forwardDepth[row - 1];
        int leftDepth = 0;
        int leftQuery = 0;
        int rightDepth = topology.maxDepth + 1;
        int rightQuery = queryLength;
        for (POASeedAnchor anchor : anchors) {
            int anchorDepth = topology.forwardDepth[anchor.graphStartRow() - 1];
            if (anchorDepth <= depth) {
                leftDepth = anchorDepth;
                leftQuery = anchor.queryStart() + 1;
            }
            if (anchorDepth >= depth) {
                rightDepth = anchorDepth;
                rightQuery = anchor.queryStart() + 1;
                break;
            }
        }
        if (rightDepth == leftDepth) return leftQuery;
        double fraction = (depth - leftDepth) / (double) (rightDepth - leftDepth);
        return (int) Math.round(leftQuery + fraction * (rightQuery - leftQuery));
    }

    private AlignmentRun runExistingDense(String query, long priorEvaluatedCells) {
        POAAlignmentPath best = poa.alignBestToGraph(query);
        long denseCells = (long) topology.nodes.size() * query.length();
        return new AlignmentRun(best.path(), best.score(), query.length(),
                best.queryInterval(), true, false,
                priorEvaluatedCells + denseCells);
    }

    private AlignmentRun run(String query, RowIntervals intervals) {
        Workspace workspace = workspaces.get();
        if (workspace == null) {
            workspace = new Workspace(intervals, topology.nodes.size());
            workspaces.set(workspace);
        } else {
            workspace.reset(intervals, topology.nodes.size());
        }
        int[] boundaryScores = workspace.boundaryScores;
        int[] boundaryPrevious = workspace.boundaryPrevious;
        for (int row = 1; row <= topology.nodes.size(); row++) {
            if (localMode || topology.predecessors[row - 1].length == 0) {
                boundaryScores[row] = 0;
                boundaryPrevious[row] = 0;
            } else {
                boundaryScores[row] = row * scores.gapPenalty();
                boundaryPrevious[row] = row - 1;
            }
        }
        for (int row = 1; row <= topology.nodes.size(); row++) {
            int[] predecessors = topology.predecessors[row - 1];
            char graphBase = topology.nodes.get(row - 1).getCH();
            for (int column = intervals.starts[row]; column <= intervals.ends[row]; column++) {
                if (intervals.ends[row] < intervals.starts[row]) continue;
                workspace.evaluatedCells++;
                int substitution = scores.score(graphBase, query.charAt(column - 1));
                PointerType substitutionType = graphBase == query.charAt(column - 1)
                        ? PointerType.MATCH : PointerType.MIS;
                int bestScore = UNREACHABLE;
                int bestRow = -1;
                int bestColumn = -1;
                PointerType bestType = PointerType.UNDEF;
                if (predecessors.length == 0) {
                    int prior = topBoundary(column - 1);
                    bestScore = prior + substitution;
                    bestRow = 0; bestColumn = column - 1; bestType = substitutionType;
                }
                for (int predecessor : predecessors) {
                    int prior = column == 1 ? boundaryScores[predecessor]
                            : workspace.score(predecessor, column - 1);
                    if (prior != UNREACHABLE && prior + substitution > bestScore) {
                        bestScore = prior + substitution;
                        bestRow = predecessor; bestColumn = column - 1;
                        bestType = substitutionType;
                    }
                }
                for (int predecessor : predecessors) {
                    int prior = workspace.score(predecessor, column);
                    if (prior != UNREACHABLE && prior + scores.gapPenalty() > bestScore) {
                        bestScore = prior + scores.gapPenalty();
                        bestRow = predecessor; bestColumn = column; bestType = PointerType.DEL;
                    }
                }
                int insertionPrior = column == 1 ? boundaryScores[row]
                        : workspace.score(row, column - 1);
                if (insertionPrior != UNREACHABLE
                        && insertionPrior + scores.gapPenalty() > bestScore) {
                    bestScore = insertionPrior + scores.gapPenalty();
                    bestRow = row; bestColumn = column - 1; bestType = PointerType.INS;
                }
                if (localMode && bestScore <= 0) {
                    bestScore = 0;
                    bestRow = 0;
                    bestColumn = 0;
                    bestType = PointerType.JUMP;
                }
                if (bestScore != UNREACHABLE) {
                    workspace.set(row, column, bestScore, bestRow, bestColumn, bestType);
                }
            }
        }
        int endpoint = -1;
        int endpointScore = UNREACHABLE;
        int endpointColumn = query.length();
        if (localMode) {
            endpoint = 0;
            endpointColumn = 0;
            endpointScore = 0;
            for (int row = 1; row <= topology.nodes.size(); row++) {
                for (int column = intervals.starts[row]; column <= intervals.ends[row]; column++) {
                    int score = workspace.score(row, column);
                    if (score > endpointScore) {
                        endpoint = row;
                        endpointColumn = column;
                        endpointScore = score;
                    }
                }
            }
        } else {
            for (int row : ((POGraphGlobal) poa).alignmentEndRowsInSelectionOrder()) {
                int score = workspace.score(row, query.length());
                if (score > endpointScore) {
                    endpoint = row;
                    endpointScore = score;
                }
            }
        }
        if (endpoint < 0) {
            return AlignmentRun.unreachable(query.length(), workspace.evaluatedCells);
        }
        List<POGNode> path = new ArrayList<POGNode>();
        boolean boundaryTouched = false;
        int row = endpoint;
        int column = endpointColumn;
        int alignmentEnd = endpointColumn;
        while (row > 0 || column > 0) {
            if (row == 0) {
                if (localMode) break;
                path.add(null); column--; continue;
            }
            if (column == 0) {
                if (localMode) break;
                row = boundaryPrevious[row]; continue;
            }
            if (!workspace.contains(row, column)) {
                return AlignmentRun.unreachable(query.length(), workspace.evaluatedCells);
            }
            if ((column == intervals.starts[row] && column > 1)
                    || (column == intervals.ends[row] && column < query.length())) {
                boundaryTouched = true;
            }
            PointerType type = workspace.type(row, column);
            if (type == PointerType.JUMP) break;
            if (type == PointerType.MATCH) path.add(topology.nodes.get(row - 1));
            else if (type == PointerType.MIS || type == PointerType.INS) path.add(null);
            int nextRow = workspace.previousRow(row, column);
            int nextColumn = workspace.previousColumn(row, column);
            if (nextRow < 0 || nextColumn < 0) {
                return AlignmentRun.unreachable(query.length(), workspace.evaluatedCells);
            }
            row = nextRow; column = nextColumn;
        }
        Collections.reverse(path);
        return new AlignmentRun(path, endpointScore, query.length(),
                new Interval(column, alignmentEnd), true,
                boundaryTouched, workspace.evaluatedCells);
    }

    private int topBoundary(int column) {
        return localMode ? 0 : column * scores.gapPenalty();
    }

    private static List<SeedToken> selectSeeds(String sequence, int k,
            SeedSelection selection, int window) {
        List<SeedToken> all = new ArrayList<SeedToken>();
        for (int position = 0; position + k <= sequence.length(); position++) {
            String kmer = sequence.substring(position, position + k).toUpperCase();
            if (isUnambiguousDna(kmer)) all.add(new SeedToken(position, kmer));
        }
        if (selection == SeedSelection.ALL_KMERS || all.size() <= 1) return all;
        List<SeedToken> selected = new ArrayList<SeedToken>();
        Set<Integer> positions = new HashSet<Integer>();
        int effectiveWindow = Math.min(window, all.size());
        for (int start = 0; start + effectiveWindow <= all.size(); start++) {
            SeedToken best = all.get(start);
            for (int index = start + 1; index < start + effectiveWindow; index++) {
                SeedToken candidate = all.get(index);
                int comparison = Integer.compare(candidate.kmer.hashCode(), best.kmer.hashCode());
                if (comparison < 0 || (comparison == 0 && candidate.position < best.position)) {
                    best = candidate;
                }
            }
            if (positions.add(best.position)) selected.add(best);
        }
        return selected;
    }

    private static boolean isUnambiguousDna(String value) {
        for (int i = 0; i < value.length(); i++) {
            char base = value.charAt(i);
            if (base != 'A' && base != 'C' && base != 'G' && base != 'T') return false;
        }
        return true;
    }

    private static long fingerprint(DefaultDirectedGraph<POGNode,POGEdge> graph) {
        long value = 1469598103934665603L;
        for (POGNode node : graph.vertexSet()) {
            value = (value ^ node.getCounter()) * 1099511628211L;
            value = (value ^ node.getCH()) * 1099511628211L;
            for (int read : node.getIDX()) value = (value ^ read) * 1099511628211L;
        }
        for (POGEdge edge : graph.edgeSet()) {
            value ^= graph.getEdgeSource(edge).getCounter() * 31L
                    + graph.getEdgeTarget(edge).getCounter();
        }
        return value;
    }

    private static final class SeedIndex {
        final Map<String,List<Occurrence>> byKmer;

        SeedIndex(Map<String,List<Occurrence>> byKmer) { this.byKmer = byKmer; }

        static SeedIndex build(Topology topology, SeededBandOptions options) {
            int maxRead = -1;
            for (POGNode node : topology.nodes) {
                for (int read : node.getIDX()) if (read > maxRead) maxRead = read;
            }
            List<List<Integer>> readRows = new ArrayList<List<Integer>>(maxRead + 1);
            for (int read = 0; read <= maxRead; read++) readRows.add(new ArrayList<Integer>());
            for (int row = 1; row <= topology.nodes.size(); row++) {
                for (int read : topology.nodes.get(row - 1).getIDX()) {
                    if (read >= 0) readRows.get(read).add(row);
                }
            }
            Map<String,List<Occurrence>> index = new HashMap<String,List<Occurrence>>();
            Map<String,Set<String>> seenPaths = new HashMap<String,Set<String>>();
            for (List<Integer> path : readRows) {
                StringBuilder sequence = new StringBuilder(path.size());
                for (int row : path) sequence.append(topology.nodes.get(row - 1).getCH());
                for (SeedToken seed : selectSeeds(sequence.toString(), options.kmerLength(),
                        options.seedSelection(), options.minimizerWindow())) {
                    int start = path.get(seed.position);
                    int end = path.get(seed.position + options.kmerLength() - 1);
                    StringBuilder signature = new StringBuilder();
                    for (int offset = 0; offset < options.kmerLength(); offset++) {
                        signature.append(path.get(seed.position + offset)).append(',');
                    }
                    Set<String> paths = seenPaths.get(seed.kmer);
                    if (paths == null) {
                        paths = new HashSet<String>(); seenPaths.put(seed.kmer, paths);
                    }
                    if (!paths.add(signature.toString())) continue;
                    List<Occurrence> occurrences = index.get(seed.kmer);
                    if (occurrences == null) {
                        occurrences = new ArrayList<Occurrence>(); index.put(seed.kmer, occurrences);
                    }
                    occurrences.add(new Occurrence(start, end));
                }
            }
            return new SeedIndex(index);
        }
    }

    private static final class Topology {
        final List<POGNode> nodes;
        final int[][] predecessors;
        final int[][] successors;
        final BitSet[] descendants;
        final int[] forwardDepth;
        final int maxDepth;

        Topology(List<POGNode> nodes, int[][] predecessors, int[][] successors,
                BitSet[] descendants, int[] forwardDepth, int maxDepth) {
            this.nodes = nodes; this.predecessors = predecessors; this.successors = successors;
            this.descendants = descendants; this.forwardDepth = forwardDepth;
            this.maxDepth = maxDepth;
        }

        boolean reaches(int sourceRow, int targetRow) {
            return sourceRow == targetRow
                    || descendants[sourceRow - 1].get(targetRow - 1);
        }

        static Topology build(DefaultDirectedGraph<POGNode,POGEdge> graph,
                POAAlignmentTopologyView view) {
            List<POGNode> nodes = new ArrayList<POGNode>(view.nodes);
            if (nodes.size() != graph.vertexSet().size()) {
                throw new IllegalStateException("POA graph is not acyclic");
            }
            Map<POGNode,Integer> rows = new HashMap<POGNode,Integer>();
            for (int i = 0; i < nodes.size(); i++) rows.put(nodes.get(i), i + 1);
            int[][] predecessors = new int[nodes.size()][];
            int[][] successors = new int[nodes.size()][];
            int[] depth = new int[nodes.size()];
            int maxDepth = 0;
            for (int i = 0; i < nodes.size(); i++) {
                POGNode node = nodes.get(i);
                predecessors[i] = view.predecessorRows[i].clone();
                successors[i] = edgeRows(graph, graph.outgoingEdgesOf(node), rows, false);
                int value = 1;
                for (int predecessor : predecessors[i]) value = Math.max(value, depth[predecessor - 1] + 1);
                depth[i] = value; maxDepth = Math.max(maxDepth, value);
            }
            BitSet[] descendants = new BitSet[nodes.size()];
            for (int i = nodes.size() - 1; i >= 0; i--) {
                descendants[i] = new BitSet(nodes.size());
                for (int successor : successors[i]) {
                    descendants[i].set(successor - 1);
                    descendants[i].or(descendants[successor - 1]);
                }
            }
            return new Topology(nodes, predecessors, successors, descendants, depth, maxDepth);
        }

        private static int[] edgeRows(DefaultDirectedGraph<POGNode,POGEdge> graph,
                Set<POGEdge> edges, Map<POGNode,Integer> rows, boolean sources) {
            int[] result = new int[edges.size()]; int index = 0;
            for (POGEdge edge : edges) result[index++] = rows.get(sources
                    ? graph.getEdgeSource(edge) : graph.getEdgeTarget(edge));
            return result;
        }
    }

    private static final class Workspace {
        private static final byte UNDEFINED_OPERATION =
                (byte) PointerType.UNDEF.ordinal();

        RowIntervals intervals;
        int[] offsets;
        int[] values;
        int[] previousRows;
        int[] previousColumns;
        byte[] operations;
        int[] boundaryScores;
        int[] boundaryPrevious;
        long evaluatedCells;

        Workspace(RowIntervals intervals, int graphRows) {
            offsets = new int[intervals.starts.length];
            values = new int[0];
            previousRows = new int[0];
            previousColumns = new int[0];
            operations = new byte[0];
            boundaryScores = new int[graphRows + 1];
            boundaryPrevious = new int[graphRows + 1];
            reset(intervals, graphRows);
        }

        void reset(RowIntervals newIntervals, int graphRows) {
            this.intervals = newIntervals;
            if (offsets.length < newIntervals.starts.length) {
                offsets = new int[newIntervals.starts.length];
            }
            if (boundaryScores.length < graphRows + 1) {
                boundaryScores = new int[graphRows + 1];
                boundaryPrevious = new int[graphRows + 1];
            }
            int size = 0;
            for (int row = 1; row < newIntervals.starts.length; row++) {
                offsets[row] = size;
                if (newIntervals.ends[row] >= newIntervals.starts[row]) {
                    size += newIntervals.ends[row] - newIntervals.starts[row] + 1;
                }
            }
            ensureCellCapacity(size);
            Arrays.fill(values, 0, size, UNREACHABLE);
            Arrays.fill(previousRows, 0, size, -1);
            Arrays.fill(previousColumns, 0, size, -1);
            Arrays.fill(operations, 0, size, UNDEFINED_OPERATION);
            evaluatedCells = 0;
        }

        private void ensureCellCapacity(int required) {
            if (values.length >= required) {
                return;
            }
            int capacity = grow(values.length, required);
            values = new int[capacity];
            previousRows = new int[capacity];
            previousColumns = new int[capacity];
            operations = new byte[capacity];
        }

        private static int grow(int current, int required) {
            int grown = current == 0 ? 1 : current;
            while (grown < required) {
                grown = grown < 1024 ? grown * 2 : grown + grown / 2;
            }
            return grown;
        }

        boolean contains(int row, int column) {
            return row > 0 && row < intervals.starts.length
                    && column >= intervals.starts[row] && column <= intervals.ends[row];
        }
        int index(int row, int column) { return offsets[row] + column - intervals.starts[row]; }
        int score(int row, int column) { return contains(row, column) ? values[index(row, column)] : UNREACHABLE; }
        void set(int row, int column, int score, int priorRow, int priorColumn, PointerType type) {
            int index = index(row, column); values[index] = score;
            previousRows[index] = priorRow; previousColumns[index] = priorColumn;
            operations[index] = (byte) type.ordinal();
        }
        int previousRow(int row, int column) { return previousRows[index(row, column)]; }
        int previousColumn(int row, int column) { return previousColumns[index(row, column)]; }
        PointerType type(int row, int column) {
            return POINTER_TYPES[operations[index(row, column)]];
        }
    }

    private static final class RowIntervals {
        final int[] starts; final int[] ends;
        RowIntervals(int[] starts, int[] ends) { this.starts = starts; this.ends = ends; }
        int activeRows() { int count = 0; for (int row = 1; row < ends.length; row++) if (ends[row] >= starts[row]) count++; return count; }
    }

    private static final class AlignmentRun {
        final List<POGNode> path; final int score; final int queryLength;
        final Interval queryInterval;
        final boolean reachable; final boolean boundaryTouched; long evaluatedCells;
        AlignmentRun(List<POGNode> path, int score, int queryLength, Interval queryInterval,
                boolean reachable, boolean boundaryTouched, long evaluatedCells) {
            this.path = path; this.score = score; this.queryLength = queryLength;
            this.queryInterval = queryInterval;
            this.reachable = reachable; this.boundaryTouched = boundaryTouched;
            this.evaluatedCells = evaluatedCells;
        }
        static AlignmentRun unreachable(int queryLength, long cells) {
            return new AlignmentRun(Collections.<POGNode>emptyList(), UNREACHABLE,
                    queryLength, new Interval(0, queryLength), false, false, cells);
        }
    }

    private static final class SeedToken {
        final int position; final String kmer;
        SeedToken(int position, String kmer) { this.position = position; this.kmer = kmer; }
    }
    private static final class Occurrence {
        final int startRow; final int endRow;
        Occurrence(int startRow, int endRow) { this.startRow = startRow; this.endRow = endRow; }
    }
    private static final class Hit {
        final int queryStart; final Occurrence occurrence; final int occurrenceCount; final int k;
        Hit(int queryStart, Occurrence occurrence, int occurrenceCount, int k) {
            this.queryStart = queryStart; this.occurrence = occurrence;
            this.occurrenceCount = occurrenceCount; this.k = k;
        }
        static final Comparator<Hit> RARITY_ORDER = new Comparator<Hit>() {
            @Override public int compare(Hit a, Hit b) {
                int c = Integer.compare(a.occurrenceCount, b.occurrenceCount);
                if (c == 0) c = Integer.compare(a.queryStart, b.queryStart);
                if (c == 0) c = Integer.compare(a.occurrence.startRow, b.occurrence.startRow);
                return c;
            }
        };
        static final Comparator<Hit> CHAIN_ORDER = new Comparator<Hit>() {
            @Override public int compare(Hit a, Hit b) {
                int c = Integer.compare(a.queryStart, b.queryStart);
                if (c == 0) c = Integer.compare(a.occurrence.startRow, b.occurrence.startRow);
                if (c == 0) c = Integer.compare(a.occurrence.endRow, b.occurrence.endRow);
                return c;
            }
        };
    }
    private static final class SeedSearch {
        final int querySeeds; final int repetitiveFiltered; final int hitCount;
        final List<POASeedAnchor> chain;
        SeedSearch(int querySeeds, int repetitiveFiltered, int hitCount, List<POASeedAnchor> chain) {
            this.querySeeds = querySeeds; this.repetitiveFiltered = repetitiveFiltered;
            this.hitCount = hitCount; this.chain = chain;
        }
    }
    private static final class IndexKey {
        final int k; final SeedSelection selection; final int window;
        IndexKey(int k, SeedSelection selection, int window) { this.k = k; this.selection = selection; this.window = window; }
        @Override public int hashCode() { return 31 * (31 * k + selection.hashCode()) + window; }
        @Override public boolean equals(Object other) {
            if (!(other instanceof IndexKey)) return false;
            IndexKey key = (IndexKey) other;
            return k == key.k && selection == key.selection && window == key.window;
        }
    }
}
