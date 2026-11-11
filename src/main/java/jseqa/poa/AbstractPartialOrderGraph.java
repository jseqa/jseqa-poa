package jseqa.poa;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.stream.Collectors;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.jgrapht.UndirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.SimpleGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import org.jgrapht.traverse.TopologicalOrderIterator;

import jseqa.alignment.AlignmentMatrixScore;
import jseqa.alignment.substitution.DnaSubstitutionMatrix;
import jseqa.alignment.substitution.SubstitutionMatrix;
import jseqa.poa.POGPtr.PointerType;
import jseqa.util.Interval;
import jseqa.util.Pair;

/**
 * Provides the skeletal implementation of Partial Order Graphs.
 *
 */
public abstract class AbstractPartialOrderGraph implements Iterable<POGNode>, PartialOrderAlignment {

	final static Logger log = LogManager.getLogger(AbstractPartialOrderGraph.class);
	private static final int COLOR_COHERENCE_SCALE = 1000000;
	
	final private Pair<Integer,Integer> originPtr = new Pair<Integer,Integer>(0,0);
	
	/**
	 * Index of sequences added to POGraph.
	 */
	protected int readIDX;

	/**
	 * Original progressive inputs, used to retain locally unaligned read ends in MSAs.
	 */
	private List<String> readSequences;
	
	/**
	 * Map from POGNode index.
	 */
	private HashMap<POGNode,Integer> nodeMap;
	
	/**
	 * List of nodes, providing a map back from index POGNode.
	 */
	protected ArrayList<POGNode> nodeList;

	/** Cached primitive topology used by repeated read-only alignments. */
	private char[] nodeBases;
	private int[][] predecessorRows;
	private int[][] predecessorColorOverlaps;
	private boolean topologyCacheValid;
	private int cachedVertexCount;
	private int cachedEdgeCount;
	private int topologyBuildCount;
	private boolean colorCoherenceTieBreaking;
	
	/**
	 * Scoring matrix, rows = graph, cols = sequence.
	 */
//	protected AlignMatrixScore alignScoreMat;
	protected AlignmentMatrixScore alignScoreMat;
	/** Cumulative color-path coherence used to resolve primary-score ties. */
	protected ReusableAlignmentMatrixScore coherenceScoreMat;
	
	/**
	 * matrix of pointers to cells.
	 */
	protected AlignPOAMatrixPtr alignPtrMat;
	
	/** matrix of backtracking pointers */
	protected AlignPOGPtrMatrix ptrMat;
	
	/** Scoring matrix to use for scoring matches/mismatches between characters */
//	protected ScoreMatrix scoreMat;
	protected SubstitutionMatrix scoreMat;
	
	/** Partial order graph */
	protected DefaultDirectedGraph<POGNode,POGEdge> g;
	
	/** Graph of aligned nodes */
	protected UndirectedGraph<POGNode,POGEdge> alignedG;
	
	/** {@link Interval} of the aligned query sequence. */
	protected Interval queryIntval;

	Map<POGEdge, Set<Integer>> edgeReadIdx;

	/** Strategy used to create nodes and attach per-read observations. */
	protected final POGNodeFactory nodeFactory;

	/**
	 * Default constructor.
	 */
	public AbstractPartialOrderGraph() {
		this(POGNodeFactory.standard());
	}

	protected AbstractPartialOrderGraph(POGNodeFactory nodeFactory) {
		this.nodeFactory = requireNodeFactory(nodeFactory);
		this.init();
		this.initMatrix();
	}
	
	/**
	 * Constructor that initializes with a single sequence.
	 * @param initseq {@code String} sequence to initialize.
	 */
	public AbstractPartialOrderGraph(String initseq) {
		this(initseq, POGNodeFactory.standard());
	}

	protected AbstractPartialOrderGraph(String initseq, POGNodeFactory nodeFactory) {
		this.nodeFactory = requireNodeFactory(nodeFactory);
		this.init();
		this.initMatrix();
		this.initGraph(initseq);
	}

	protected AbstractPartialOrderGraph(Collection<String> alignment) {
		this(alignment, POGNodeFactory.standard());
	}

	protected AbstractPartialOrderGraph(Collection<String> alignment,
			POGNodeFactory nodeFactory) {
		this.nodeFactory = requireNodeFactory(nodeFactory);
		this.init();
		this.initMatrix();
		this.initGraphFromAlignment(alignment);
	}
	
	/**
	 * Constructor that specifies a particular scoring matrix to use
	 * @param sm {@link SubstitutionMatrix} to use for scoring the alignments
	 */
	public AbstractPartialOrderGraph(SubstitutionMatrix sm) {
		this(sm, POGNodeFactory.standard());
	}

	protected AbstractPartialOrderGraph(SubstitutionMatrix sm, POGNodeFactory nodeFactory) {
		this.nodeFactory = requireNodeFactory(nodeFactory);
		this.init(sm);
		this.initMatrix();
	}
	
	/**
	 * Constructor that initializes the graph with the sequence and uses the
	 * provided {@code ScoreMatrix} to score alignments.
	 * @param initSeq {@code String} sequence to initialize the graph
	 * @param sm {@code ScoreMatrix} to score alignments
	 */
	public AbstractPartialOrderGraph(String initSeq, SubstitutionMatrix sm) {
		this(initSeq, sm, POGNodeFactory.standard());
	}

	protected AbstractPartialOrderGraph(String initSeq, SubstitutionMatrix sm,
			POGNodeFactory nodeFactory) {
		this.nodeFactory = requireNodeFactory(nodeFactory);
		this.init(sm);
		this.initMatrix();
		this.initGraph(initSeq);
	}

	protected AbstractPartialOrderGraph(Collection<String> alignment,
			SubstitutionMatrix sm) {
		this(alignment, sm, POGNodeFactory.standard());
	}

	protected AbstractPartialOrderGraph(Collection<String> alignment,
			SubstitutionMatrix sm, POGNodeFactory nodeFactory) {
		this.nodeFactory = requireNodeFactory(nodeFactory);
		this.init(sm);
		this.initMatrix();
		this.initGraphFromAlignment(alignment);
	}
	
	/**
	 * Copy constructor that performs an independent copy of graph data structures.
	 * @param otherPOG
	 */
	public AbstractPartialOrderGraph(AbstractPartialOrderGraph otherPOG) {
		this.nodeFactory = otherPOG.nodeFactory;
		this.readSequences = new ArrayList<String>(otherPOG.readSequences);
		this.colorCoherenceTieBreaking = otherPOG.colorCoherenceTieBreaking;
		alignScoreMat = new ReusableAlignmentMatrixScore(1, 1);
		coherenceScoreMat = new ReusableAlignmentMatrixScore(1, 1);
		alignPtrMat = new AlignPOAMatrixPtr(1, 1);
		this.ptrMat = new AlignPOGPtrMatrix(1, 1);

		this.g = new DefaultDirectedGraph<POGNode,POGEdge>(POGEdge.class);
		this.alignedG = new SimpleGraph<POGNode,POGEdge>(POGEdge.class);
		Map<POGNode,POGNode> copiedNodes = new HashMap<>();
		for (POGNode node : otherPOG.g.vertexSet()) {
			POGNode copy = this.nodeFactory.copy(node);
			copiedNodes.put(node, copy);
			this.g.addVertex(copy);
		}
		for (POGEdge edge : otherPOG.g.edgeSet()) {
			POGEdge copy = new POGEdge(edge.getEdgeType());
			this.g.addEdge(copiedNodes.get(otherPOG.g.getEdgeSource(edge)),
					copiedNodes.get(otherPOG.g.getEdgeTarget(edge)), copy);
		}
		for (POGNode node : otherPOG.alignedG.vertexSet()) {
			this.alignedG.addVertex(copiedNodes.get(node));
		}
		for (POGEdge edge : otherPOG.alignedG.edgeSet()) {
			this.alignedG.addEdge(copiedNodes.get(otherPOG.alignedG.getEdgeSource(edge)),
					copiedNodes.get(otherPOG.alignedG.getEdgeTarget(edge)),
					new POGEdge(edge.getEdgeType()));
		}

		scoreMat = otherPOG.scoreMat;
		readIDX = otherPOG.readIDX;
		this.queryIntval = otherPOG.queryIntval == null ? null
				: new Interval(otherPOG.queryIntval.start(), otherPOG.queryIntval.end());
		this.edgeReadIdx = new HashMap<>();
		for (POGEdge edge : otherPOG.g.edgeSet()) {
			POGEdge copiedEdge = this.g.getEdge(copiedNodes.get(otherPOG.g.getEdgeSource(edge)),
					copiedNodes.get(otherPOG.g.getEdgeTarget(edge)));
			Set<Integer> indices = otherPOG.edgeReadIdx.get(edge);
			this.edgeReadIdx.put(copiedEdge,
					indices == null ? new HashSet<Integer>() : new HashSet<Integer>(indices));
		}
		this.computeInternalDS();
		this.initMatrix();
	}

	private static POGNodeFactory requireNodeFactory(POGNodeFactory factory) {
		if (factory == null) {
			throw new IllegalArgumentException("node factory must not be null");
		}
		return factory;
	}
	
	/**
	 * Adds a sequence to the partial order alignment graph
	 * @param seq {@code String} sequence to add to the graph
	 */
	public void add(String seq) {
		validateSequence(seq);
		if(this.g.vertexSet().size() == 0) {
			this.initGraph(seq);
		} else {
			this.alignSequence(seq);
			this.backtrack(seq);
			readIDX++;
			readSequences.add(seq);
		}
	}
	
	/**
	 * Adds all sequences in the collection to the partial order alignment graph. Is merely multiple 
	 * calls to {@link #add(String)}.
	 * @param seqs {@code Collection} of {@code String} sequences to add
	 */
	public void addAll(Collection<String> seqs) {		
		for(String seq : seqs) {
			if(this.g.vertexSet().size() == 0) {
				this.initGraph(seq);
			} else {
				this.add(seq);				
			}
		}
	}

	public Set<Integer> getEdgeIdx(POGEdge e) {
		return this.edgeReadIdx.get(e);
	}

	/**
	 * Enables a cumulative color-coherence objective for resolving paths tied on
	 * their primary alignment score. The secondary score sums the fraction of
	 * each node's reference colors retained across the selected incoming edge.
	 *
	 * @param enabled whether to prefer color-coherent tied paths
	 */
	public void setColorCoherenceTieBreaking(boolean enabled) {
		this.colorCoherenceTieBreaking = enabled;
	}

	/**
	 * Align the sequence to the POGraph without adding nodes/edges, and return the backtracked list of 
	 * {@link POGNode} instances representing the path of the aligned sequence.
	 * @param seq {@code String} sequence to align to the POGraph
	 * @return List of {@code POGNodes} representing the path of the aligned sequence
	 */
	public List<POGNode> alignToGraph(String seq) {
		validateSequence(seq);
		this.alignSequence(seq);
		return this.backtrackPath(seq);
	}

	@Override
	public POAAlignmentPath alignBestToGraph(String seq) {
		validateSequence(seq);
		this.alignSequence(seq);
		Pair<Integer,Integer> endCell = this.backtrackFromCell(seq.length());
		int score = this.alignScoreMat.get(endCell.first(), endCell.second());
		List<POGNode> path = this.backtrackPath(seq, endCell);
		return new POAAlignmentPath(path, this.queryIntval, score);
	}

	@Override
	public List<POAAlignmentPath> alignToGraph(String seq, int k) {
		validateSequence(seq);
		if (k < 1) {
			throw new IllegalArgumentException("k must be at least one");
		}
		computeInternalDS();
		prepareMatrices(seq.length());
		KBestAlignmentMatrix matrix = new KBestAlignmentMatrix(
				nodeList.size() + 1, seq.length() + 1, k);
		seedKBestBoundaries(matrix);
		for (int row = 1; row <= nodeList.size(); row++) {
			int[] predecessors = predecessorRows[row - 1];
			char graphBase = nodeBases[row - 1];
			for (int column = 1; column <= seq.length(); column++) {
				fillKBestCell(matrix, row, column, seq.charAt(column - 1),
						graphBase, predecessors);
			}
		}
		List<KBestEndpoint> endpoints = selectKBestEndpoints(matrix, seq.length(), k);
		List<POAAlignmentPath> paths = new ArrayList<POAAlignmentPath>(endpoints.size());
		Set<String> seen = new HashSet<String>();
		for (KBestEndpoint endpoint : endpoints) {
			POAAlignmentPath path = backtrackKBestPath(matrix, endpoint);
			String signature = pathSignature(path);
			if (seen.add(signature)) {
				paths.add(path);
			}
		}
		if (!paths.isEmpty()) {
			this.queryIntval = paths.get(0).queryInterval();
		}
		return paths;
	}

	private void seedKBestBoundaries(KBestAlignmentMatrix matrix) {
		for (int row = 0; row < matrix.rows(); row++) {
			seedKBestBoundaryCell(matrix, row, 0);
		}
		for (int column = 1; column < matrix.columns(); column++) {
			seedKBestBoundaryCell(matrix, 0, column);
		}
	}

	private void seedKBestBoundaryCell(KBestAlignmentMatrix matrix, int row, int column) {
		POGPtr pointer = ptrMat.get(row, column);
		PointerType type = row == 0 && column == 0 ? PointerType.UNDEF : pointer.type();
		matrix.set(row, column, 0, alignScoreMat.get(row, column), 0,
				pointer.ptr().first(), pointer.ptr().second(), 0, type);
	}

	private void fillKBestCell(KBestAlignmentMatrix matrix, int row, int column,
			char sequenceBase, char graphBase, int[] predecessors) {
		KBestCellAccumulator best = new KBestCellAccumulator(matrix.ranks());
		// Local alignment replaces every non-positive cell with a restart. Adding
		// it first preserves that tie order while positive candidates still win.
		if (allowsKBestLocalRestart()) {
			best.add(0, 0, 0, 0, 0, PointerType.JUMP);
		}
		int substitutionScore = scoreMat.score(graphBase, sequenceBase);
		PointerType substitutionType = graphBase == sequenceBase
				? PointerType.MATCH : PointerType.MIS;
		if (predecessors.length == 0) {
			addKBestCandidates(best, matrix, 0, column - 1,
					substitutionScore, 0, substitutionType);
		}
		for (int predecessor : predecessors) {
			int coherence = colorCoherenceTieBreaking
					? predecessorColorCoherence(row, predecessor) : 0;
			addKBestCandidates(best, matrix, predecessor, column - 1,
					substitutionScore, coherence, substitutionType);
		}
		// Keep substitutions ahead of deletions, matching fillCell's legacy
		// strict-greater tie behavior for the rank-one state.
		for (int predecessor : predecessors) {
			int coherence = colorCoherenceTieBreaking
					? predecessorColorCoherence(row, predecessor) : 0;
			addKBestCandidates(best, matrix, predecessor, column,
					scoreMat.gapPenalty(), coherence, PointerType.DEL);
		}
		addKBestCandidates(best, matrix, row, column - 1,
				scoreMat.gapPenalty(), 0, PointerType.INS);
		best.writeTo(matrix, row, column);
	}

	private void addKBestCandidates(KBestCellAccumulator best,
			KBestAlignmentMatrix matrix, int previousRow,
			int previousColumn, int scoreDelta, int coherenceDelta,
			PointerType operation) {
		for (int rank = 0; rank < matrix.ranks(); rank++) {
			int previousScore = matrix.score(previousRow, previousColumn, rank);
			if (previousScore == KBestAlignmentMatrix.UNREACHABLE) {
				continue;
			}
			best.add(previousScore + scoreDelta,
					matrix.coherence(previousRow, previousColumn, rank) + coherenceDelta,
					previousRow, previousColumn, rank, operation);
		}
	}

	/** Whether an interior cell can restart a local alignment at zero. */
	protected boolean allowsKBestLocalRestart() {
		return false;
	}

	/** Whether a cell can terminate a k-best alignment for this alignment mode. */
	protected abstract boolean isKBestEndCell(int row, int column, int queryLength);

	private List<KBestEndpoint> selectKBestEndpoints(KBestAlignmentMatrix matrix,
			int queryLength, int k) {
		KBestEndpointAccumulator best = new KBestEndpointAccumulator(k);
		for (int row = 0; row < matrix.rows(); row++) {
			for (int column = 0; column < matrix.columns(); column++) {
				if (!isKBestEndCell(row, column, queryLength)) {
					continue;
				}
				for (int rank = 0; rank < matrix.ranks(); rank++) {
					int score = matrix.score(row, column, rank);
					if (score != KBestAlignmentMatrix.UNREACHABLE) {
						best.add(new KBestEndpoint(row, column, rank, score,
								matrix.coherence(row, column, rank)));
					}
				}
			}
		}
		return best.results();
	}

	private POAAlignmentPath backtrackKBestPath(
			KBestAlignmentMatrix matrix, KBestEndpoint endpoint) {
		LinkedList<POGNode> path = new LinkedList<POGNode>();
		int row = endpoint.row;
		int column = endpoint.column;
		int rank = endpoint.rank;
		while (row != 0 || column != 0) {
			PointerType operation = matrix.operation(row, column, rank);
			if (operation == PointerType.JUMP || operation == PointerType.UNDEF) {
				break;
			}
			if (operation == PointerType.MATCH) {
				path.addFirst(nodeList.get(row - 1));
			} else if (operation == PointerType.MIS || operation == PointerType.INS) {
				path.addFirst(null);
			}
			int previousRow = matrix.previousRow(row, column, rank);
			int previousColumn = matrix.previousColumn(row, column, rank);
			int previousRank = matrix.previousRank(row, column, rank);
			if (previousRow < 0 || previousColumn < 0 || previousRank < 0) {
				break;
			}
			row = previousRow;
			column = previousColumn;
			rank = previousRank;
		}
		return new POAAlignmentPath(path, new Interval(column, endpoint.column), endpoint.score);
	}

	private static String pathSignature(POAAlignmentPath path) {
		StringBuilder signature = new StringBuilder();
		signature.append(path.queryInterval().start()).append(':')
				.append(path.queryInterval().end()).append(':');
		for (POGNode node : path.path()) {
			signature.append(node == null ? "-" : Integer.toString(node.getCounter())).append(',');
		}
		return signature.toString();
	}

	private static final class KBestCellAccumulator {
		private final int[] scores;
		private final int[] coherence;
		private final int[] previousRows;
		private final int[] previousColumns;
		private final int[] previousRanks;
		private final PointerType[] operations;

		private KBestCellAccumulator(int k) {
			scores = new int[k];
			coherence = new int[k];
			previousRows = new int[k];
			previousColumns = new int[k];
			previousRanks = new int[k];
			operations = new PointerType[k];
			java.util.Arrays.fill(scores, KBestAlignmentMatrix.UNREACHABLE);
			java.util.Arrays.fill(previousRows, -1);
			java.util.Arrays.fill(previousColumns, -1);
			java.util.Arrays.fill(previousRanks, -1);
			java.util.Arrays.fill(operations, PointerType.UNDEF);
		}

		private void add(int score, int colorCoherence, int previousRow,
				int previousColumn, int previousRank, PointerType operation) {
			for (int index = 0; index < scores.length; index++) {
				if (scores[index] == score && coherence[index] == colorCoherence
						&& previousRows[index] == previousRow
						&& previousColumns[index] == previousColumn
						&& previousRanks[index] == previousRank
						&& operations[index] == operation) {
					return;
				}
				if (score > scores[index]
						|| (score == scores[index] && colorCoherence > coherence[index])) {
					for (int shift = scores.length - 1; shift > index; shift--) {
						copy(shift - 1, shift);
					}
					scores[index] = score;
					coherence[index] = colorCoherence;
					previousRows[index] = previousRow;
					previousColumns[index] = previousColumn;
					previousRanks[index] = previousRank;
					operations[index] = operation;
					return;
				}
			}
		}

		private void copy(int source, int target) {
			scores[target] = scores[source];
			coherence[target] = coherence[source];
			previousRows[target] = previousRows[source];
			previousColumns[target] = previousColumns[source];
			previousRanks[target] = previousRanks[source];
			operations[target] = operations[source];
		}

		private void writeTo(KBestAlignmentMatrix matrix, int row, int column) {
			for (int rank = 0; rank < scores.length; rank++) {
				if (scores[rank] != KBestAlignmentMatrix.UNREACHABLE) {
					matrix.set(row, column, rank, scores[rank], coherence[rank],
							previousRows[rank], previousColumns[rank], previousRanks[rank],
							operations[rank]);
				}
			}
		}
	}

	private static final class KBestEndpoint {
		private final int row;
		private final int column;
		private final int rank;
		private final int score;
		private final int coherence;

		private KBestEndpoint(int row, int column, int rank, int score, int coherence) {
			this.row = row;
			this.column = column;
			this.rank = rank;
			this.score = score;
			this.coherence = coherence;
		}
	}

	private static final class KBestEndpointAccumulator {
		private final KBestEndpoint[] endpoints;

		private KBestEndpointAccumulator(int k) {
			endpoints = new KBestEndpoint[k];
		}

		private void add(KBestEndpoint candidate) {
			for (int index = 0; index < endpoints.length; index++) {
				KBestEndpoint current = endpoints[index];
				if (current != null && current.row == candidate.row
						&& current.column == candidate.column && current.rank == candidate.rank) {
					return;
				}
				if (current == null || candidate.score > current.score
						|| (candidate.score == current.score
								&& candidate.coherence > current.coherence)) {
					for (int shift = endpoints.length - 1; shift > index; shift--) {
						endpoints[shift] = endpoints[shift - 1];
					}
					endpoints[index] = candidate;
					return;
				}
			}
		}

		private List<KBestEndpoint> results() {
			List<KBestEndpoint> results = new ArrayList<KBestEndpoint>(endpoints.length);
			for (KBestEndpoint endpoint : endpoints) {
				if (endpoint != null) {
					results.add(endpoint);
				}
			}
			return results;
		}
	}
	
	/**
	 * Compute the multiple sequence alignment (MSA) from the partial order graph.
	 * @return {@code List<String>} of each row of the MSA
	 */
	public List<String> getAlignment() {
		List<List<POGNode>> columns = this.computeAlignmentColumns();
		int numCols = columns.size();
//		log.debug("NODES:\t" + numNodes + "\tCOLS:\t" + numCols);
		
		int numReads = this.readIDX;
		List<char[]> algn = new ArrayList<char[]>(numReads);
		for(int i = 0; i < numReads; i++) {
			char[] row = new char[numCols];
			java.util.Arrays.fill(row, '-');
			algn.add(row);
		}

		for (int columnIndex = 0; columnIndex < columns.size(); columnIndex++) {
			boolean[] occupiedReads = new boolean[numReads];
			for(POGNode alignedNode : columns.get(columnIndex)) {
				char nodeCH = alignedNode.getCH();
				List<Integer> nodeIDX = alignedNode.getIDX();
				for(Integer idx : nodeIDX) {
					if(idx < 0) { continue; }
					if (idx >= numReads) {
						throw new IllegalStateException("node contains invalid read index " + idx);
					}
					if (occupiedReads[idx]) {
						throw new IllegalStateException("read " + idx
								+ " occurs more than once in alignment column " + columnIndex);
					}
					occupiedReads[idx] = true;
					algn.get(idx)[columnIndex] = nodeCH;
				}
			}
		}
		
		ArrayList<String> alignLst = new ArrayList<String>(numReads);
		for(int i = 0; i < numReads; i++) {
			alignLst.add(i, new String(algn.get(i)));
		}
		
		return this.addUnalignedReadEnds(alignLst);
	}

	/** Add progressive local-alignment prefixes/suffixes as read-specific MSA columns. */
	private List<String> addUnalignedReadEnds(List<String> graphRows) {
		if (graphRows.size() != readSequences.size()) {
			throw new IllegalStateException("read sequence metadata does not match graph rows");
		}
		List<StringBuilder> rows = new ArrayList<StringBuilder>(graphRows.size());
		for (String row : graphRows) {
			rows.add(new StringBuilder(row));
		}
		for (int read = 0; read < rows.size(); read++) {
			String graphSequence = rows.get(read).toString().replace("-", "");
			String original = readSequences.get(read);
			if (graphSequence.equals(original)) {
				continue;
			}
			if (!allowsUnalignedQueryEnds()) {
				throw new IllegalStateException("graph path for read " + read
						+ " does not preserve its complete input sequence");
			}
			int alignedStart = graphSequence.isEmpty() ? 0 : original.indexOf(graphSequence);
			if (alignedStart < 0) {
				throw new IllegalStateException("graph path for read " + read
						+ " is not an ordered subsequence of its input");
			}
			String prefix = original.substring(0, alignedStart);
			String suffix = original.substring(alignedStart + graphSequence.length());
			int firstResidue = firstResidueColumn(rows.get(read));
			insertReadOnlyColumns(rows, read, firstResidue, prefix);
			int lastResidue = lastResidueColumn(rows.get(read));
			insertReadOnlyColumns(rows, read, lastResidue + 1, suffix);
		}
		List<String> completeRows = new ArrayList<String>(rows.size());
		for (StringBuilder row : rows) {
			completeRows.add(row.toString());
		}
		return completeRows;
	}

	/** Whether this alignment mode can legitimately omit query prefixes or suffixes. */
	protected boolean allowsUnalignedQueryEnds() {
		return false;
	}

	private static void insertReadOnlyColumns(List<StringBuilder> rows, int read,
			int column, String residues) {
		for (int offset = 0; offset < residues.length(); offset++) {
			int insertionColumn = column + offset;
			for (StringBuilder row : rows) {
				row.insert(insertionColumn, '-');
			}
			rows.get(read).setCharAt(insertionColumn, residues.charAt(offset));
		}
	}

	private static int firstResidueColumn(StringBuilder row) {
		for (int column = 0; column < row.length(); column++) {
			if (row.charAt(column) != '-') {
				return column;
			}
		}
		return row.length();
	}

	private static int lastResidueColumn(StringBuilder row) {
		for (int column = row.length() - 1; column >= 0; column--) {
			if (row.charAt(column) != '-') {
				return column;
			}
		}
		return -1;
	}

	public DefaultDirectedGraph<POGNode,POGEdge> getGraph() {
		// The returned graph is mutable, so rebuild cached topology before reuse.
		this.invalidateTopologyCache();
		return this.g;
	}

	/**
	 * Compute the consensus sequence, as represented by the heaviest path
	 * @return {@code String} consensus sequence of the heaviest path.
	 */
	public String consensus() {
		if (g.vertexSet().isEmpty()) {
			return "";
		}
		
		TopologicalOrderIterator<POGNode,POGEdge> titer = 
				new TopologicalOrderIterator<POGNode,POGEdge>(g);
		HashMap<POGNode,Integer> nodeScores = new HashMap<POGNode,Integer>(g.vertexSet().size());
		HashMap<POGNode,POGNode> nodePtrs = new HashMap<POGNode,POGNode>(g.vertexSet().size());
//		POGNode initNode = titer.next();
		while(titer.hasNext()) {
			POGNode currN = titer.next();
			// get nodes of incoming edges 
			Set<POGEdge> adjE = g.incomingEdgesOf(currN);
			// get all the previous nodes
			List<POGNode> prevNodes = adjE.stream()
										.map(e -> g.getEdgeSource(e))
										.collect(Collectors.toList());
//			log.debug("CONSENSUS\t" + currN + "\t" + prevNodes + "\t" + nodeScores);
			// set score as max{incoming}+#reads_in_currNode
			if(prevNodes.size() == 0) {
				nodeScores.put(currN, currN.getNumReads());
				nodePtrs.put(currN, null);
			} else {				
				// find the previous node with most reads
				POGNode bestPrev = prevNodes.stream()
											//.max(Comparator.comparing(n -> ((POGNode)n).getNumReads()))
											.max(Comparator.comparing(n -> nodeScores.get(n)))
											.get();
				int currScore = nodeScores.get(bestPrev) + currN.getNumReads();
				// assign score to current node
				nodeScores.put(currN, currScore);
				nodePtrs.put(currN, bestPrev);
			}
		}
		// find node with highest score
		POGNode maxEndNode = nodeScores.keySet()
									.stream()
									.max(Comparator.comparing(n -> nodeScores.get(n)))
									.get();
			
		// backtrack from maximal scoring node
		StringBuilder consensusSeq = new StringBuilder();
		POGNode prevNode = maxEndNode;
		while(prevNode != null) {
			consensusSeq.append(prevNode.getCH());
//			log.debug("BACKTRACK:\t" + prevNode + "\t" + consensusSeq);
			prevNode = nodePtrs.get(prevNode);
		}
		
		return consensusSeq.reverse().toString();
	}
	
//	/**
//	 * Get the path thru the graph for a specific read index
//	 * @param readIdx
//	 * @return
//	 */
//	public List<POGNode> getPathOfReadIdx(int readIdx) {
//		TopologicalOrderIterator<POGNode,POGEdge> titer = 
//				new TopologicalOrderIterator<POGNode,POGEdge>(g);
//		POGNode startNode = null;
//		int i = 1;		
//		while(titer.hasNext()) {
//			POGNode currN = titer.next();
//			if(startNode == null && currN.getIDX().contains(readIdx)) {
//				startNode = currN;
//			}
//		}
//		List<POGNode> readPath = new ArrayList<POGNode>();
//		readPath.add(startNode);
//	
////		POGNode currNode = startNode;
////		while(currNode != null && currNode.getIDX().contains(readIdx)) {
////			Set<POGEdge> outEdges = this.g.outgoingEdgesOf(currNode);
////			List<POGNode> adjNodes = outEdges.stream().map(e -> this.g.getEdgeTarget(e)).collect(Collectors.toList());
////			Optional<POGNode> adjNode = adjNodes.stream().filter(n -> n.getIDX().contains(readIdx)).findFirst();
////			if(!adjNode.isPresent()) {
////				return readPath;
////			}
////			readPath.add(adjNode.get());
////			currNode = adjNode.get();
////		}
//
//		
//		return readPath;
//	}
	
	/**
	 * Output graph to a DOT file
	 * @param dotFilename {@code String} DOT filename
	 */
	public void toDOT(String dotFilename) {
		this.toDOT(dotFilename, false);
	}

	/**
	 * Output graph to a DOT file
	 * @param dotFilename {@code String} DOT filename
	 * @param showNodeLabels boolean to show nodes with labels or not.
	 */
	public void toDOT(String dotFilename, boolean showNodeLabels) {
		this.computeInternalDS();
		try(PrintWriter pout = new PrintWriter(dotFilename)) {
			pout.println("digraph g {");
			for(int i = 0; i < this.nodeList.size(); i++) {
//				pout.println( i + " [label = \"" + i + ": "+ this.nodeList.get(i).toString() + "\"];");
				double nodeSize = (1.0+Math.ceil(Math.log((double)this.nodeList.get(i).getNumReads()))) / 2.0;
				String nodeStr = showNodeLabels
						? (" [shape=circle, width = " + nodeSize + ",\t label=\"" + this.nodeList.get(i).getCH() + "\", fontsize=\"18pt\"];")
						: " [shape=point, width = " + nodeSize + "];";
				pout.println(i + nodeStr);
				//pout.println( i + " [shape=point, width = " + nodeSize + "];");
//				pout.println( i + " [shape=circle, width = " + nodeSize + "\t, label=\"" + this.nodeList.get(i).getCH() + "\"];");
			}
			for(POGEdge edge : this.g.edgeSet()) {
				POGNode src = this.g.getEdgeSource(edge);
				POGNode dest = this.g.getEdgeTarget(edge);
				if(this.nodeMap.containsKey(src) && this.nodeMap.containsKey(dest)) {
					int src_idx = this.nodeMap.get(src);
					int dest_idx = this.nodeMap.get(dest);
					pout.println( src_idx + " -> " + dest_idx + ";");
				}
			}
			pout.println("}");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Output the multiple sequence alignment
	 */
	public String toString() {
		String str = "";
		for(String row : this.getAlignment()) {
			str += row + "\n";			
		}
		return str;
	}
	
	/**
	 * Initializes the score matrix and backtrack pointer matrix
	 */
	abstract protected void initMatrix();
	
	/**
	 * Initialize score and backtrack pointer matrices, along with the POGraph. Uses
	 * {@link DnaSubstitutionMatrix} as the default.
	 */
	protected void init() {
		this.init(new DnaSubstitutionMatrix(-4, -4));
	}
	
	/**
	 * Initialize score and backtrack pointer matrices, along with the POGraph
	 * @param sm {@code ScoreMatrix} to use for scoring
	 */
	protected void init(SubstitutionMatrix sm) {
		g = new DefaultDirectedGraph<POGNode,POGEdge>(POGEdge.class);
		alignedG = new SimpleGraph<POGNode,POGEdge>(POGEdge.class);
		alignScoreMat = new ReusableAlignmentMatrixScore(1, 1);
		coherenceScoreMat = new ReusableAlignmentMatrixScore(1, 1);
		alignPtrMat = new AlignPOAMatrixPtr(1, 1);
		this.ptrMat = new AlignPOGPtrMatrix(1, 1);
		scoreMat = sm;
		readIDX = 0;
		readSequences = new ArrayList<String>();
		this.edgeReadIdx = new HashMap<>();
	}
	
	/**
	 * initialize the graph with a single sequence, making a chain
	 * @param seq {@code String} initial sequence
	 */
	protected void initGraph(String seq) {
		validateSequence(seq);
		this.invalidateTopologyCache();
		POGNode prevNode = null;
		for(int i = 0; i < seq.length(); i++) {
			POGNode currNode = this.nodeFactory.create(seq.charAt(i), readIDX, i);
			g.addVertex( currNode );
			if(prevNode != null) {
				POGEdge e = g.addEdge(prevNode, currNode);
				this.safeEdgeAddIdx(e, readIDX);
			}
			alignedG.addVertex(currNode);
			prevNode = currNode;
		}	
		readIDX++;
		readSequences.add(seq);
	}	

	/** Build a graph that preserves the columns of an existing MSA exactly. */
	protected void initGraphFromAlignment(Collection<String> alignment) {
		List<String> rows = validateAlignment(alignment);
		this.invalidateTopologyCache();
		POGNode[] previousByRead = new POGNode[rows.size()];
		int[] positionByRead = new int[rows.size()];
		int width = rows.get(0).length();

		for (int column = 0; column < width; column++) {
			Map<Character,POGNode> nodesByResidue = new LinkedHashMap<Character,POGNode>();
			for (int read = 0; read < rows.size(); read++) {
				char residue = rows.get(read).charAt(column);
				if (residue == '-') {
					continue;
				}

				POGNode node = nodesByResidue.get(residue);
				if (node == null) {
					node = this.nodeFactory.create(residue, read, positionByRead[read]);
					nodesByResidue.put(residue, node);
					g.addVertex(node);
					alignedG.addVertex(node);
					for (POGNode alternative : nodesByResidue.values()) {
						if (alternative != node) {
							alignedG.addEdge(alternative, node);
						}
					}
				} else {
					this.nodeFactory.addReadObservation(node, read, positionByRead[read]);
				}

				POGNode previous = previousByRead[read];
				if (previous != null) {
					POGEdge edge = g.getEdge(previous, node);
					if (edge == null) {
						edge = g.addEdge(previous, node);
					}
					this.safeEdgeAddIdx(edge, read);
				}
				previousByRead[read] = node;
				positionByRead[read]++;
			}
		}

		this.readIDX = rows.size();
		for (String row : rows) {
			readSequences.add(row.replace("-", ""));
		}
		this.computeInternalDS();
	}

	private static List<String> validateAlignment(Collection<String> alignment) {
		if (alignment == null || alignment.isEmpty()) {
			throw new IllegalArgumentException("alignment must contain at least one row");
		}
		List<String> rows = new ArrayList<String>(alignment);
		if (rows.get(0) == null || rows.get(0).isEmpty()) {
			throw new IllegalArgumentException("alignment rows must not be null or empty");
		}
		int width = rows.get(0).length();
		for (String row : rows) {
			if (row == null || row.length() != width) {
				throw new IllegalArgumentException("all alignment rows must have equal width");
			}
			if (row.chars().allMatch(ch -> ch == '-')) {
				throw new IllegalArgumentException("alignment rows must contain a residue");
			}
		}
		for (int column = 0; column < width; column++) {
			boolean containsResidue = false;
			for (String row : rows) {
				if (row.charAt(column) != '-') {
					containsResidue = true;
					break;
				}
			}
			if (!containsResidue) {
				throw new IllegalArgumentException("alignment contains an all-gap column at " + column);
			}
		}
		return rows;
	}

	protected void safeEdgeAddIdx(POGEdge e, int idx) {
		if (e == null) {
			return;
		}
		if (!this.edgeReadIdx.containsKey(e)) {
			this.edgeReadIdx.put(e, new HashSet<>());
		}
		this.edgeReadIdx.get(e).add(idx);
	}

	private static void validateSequence(String sequence) {
		if (sequence == null) {
			throw new IllegalArgumentException("sequence must not be null");
		}
		if (sequence.isEmpty()) {
			throw new IllegalArgumentException("sequence must not be empty");
		}
	}

	/**
	 * Return the cell in the score matrix to start backtracking from, with first => node IDX, and second => seq IDX.
	 * @param maxSeqLen 
	 * @return
	 */
	abstract protected Pair<Integer,Integer> backtrackFromCell(int maxSeqLen);

	/** Compare completed DP cells by primary score, then cumulative color coherence. */
	protected boolean isBetterAlignmentCell(int row, int column,
			int bestRow, int bestColumn) {
		if (bestRow < 0 || bestColumn < 0) {
			return true;
		}
		int score = alignScoreMat.get(row, column);
		int bestScore = alignScoreMat.get(bestRow, bestColumn);
		return score > bestScore || (score == bestScore
				&& coherenceScoreMat.get(row, column)
				> coherenceScoreMat.get(bestRow, bestColumn));
	}
	
	abstract protected POGPtr backtrackFromPtr(int maxSeqLen);
	
	/**
	 * 
	 * @param seq
	 * @param prevNode
	 * @param currPos
	 * @param s1
	 * @param s2
	 */
	protected void backtrackBuild(String seq,
									POGNode prevNode,
									Pair<Integer,Integer> currPos, 
									String s1,
									String s2) {
		
		POGPtr currPtr = this.ptrMat.get(currPos.first(), currPos.second());
//		if (currPtr.ptr().equals(originPtr)) {
//		if (currPos.equals(originPtr)) {
		if(currPos.equals(originPtr) || currPtr.type() == PointerType.JUMP) {
			return;
		}
		
//		POGPtr currPtr = this.ptrMat.get(currPos.first(), currPos.second());
		
		int i = currPos.first();
		int j = currPos.second();

		boolean is_del = currPtr.type() == PointerType.DEL;
		boolean is_ins = currPtr.type() == PointerType.INS;
		
		POGNode currNode = i > 0 ? nodeList.get(i-1) : null;
		char nodeCH = (i == 0 ? '-' : nodeList.get(i-1).getCH());
		char seqCH = (j == 0 ? '-' : seq.charAt(j-1));
		
		s1 += (is_ins ? "-" : nodeCH);
		s2 += (is_del ? "-" : seqCH);

		// if match, add index to node
		if(currPtr.type() == PointerType.MATCH) {
			this.nodeFactory.addReadObservation(currNode, readIDX, j - 1);
//			log.trace("MAKE EDGE?:\t" + (prevNode != null) +  !g.containsEdge(currNode, prevNode));
			if (prevNode != null && !g.containsEdge(currNode, prevNode)) {
				g.addEdge(currNode, prevNode);
			}
			if (prevNode != null) {
				POGEdge e = g.getEdge(currNode, prevNode);
				this.safeEdgeAddIdx(e, readIDX);
			}

			prevNode = currNode;
		} else if(currPtr.type() == PointerType.DEL) {
			// don't add vertex, keep prevNode the same
		} else if(currPtr.type() == PointerType.INS) {
			POGNode newNode = this.nodeFactory.create(seqCH, this.readIDX, j - 1);
			g.addVertex(newNode);
			if(prevNode != null) {
				g.addEdge(newNode, prevNode);
			}
			if (prevNode != null) {
				POGEdge e = g.getEdge(newNode, prevNode);
				this.safeEdgeAddIdx(e, readIDX);
			}
			if(log.isTraceEnabled()) {
				log.trace("INSERTION:\t" + newNode + " <- " + prevNode);
			}
			prevNode = newNode;
		}
		// if mismatch, add node and connect it to all prev nodes
		else {
			POGNode newNode = this.nodeFactory.create(seqCH, this.readIDX, j - 1);
			g.addVertex(newNode);
			if(prevNode != null) {
				POGEdge e = g.addEdge(newNode, prevNode);
				this.safeEdgeAddIdx(e, readIDX);
				if(log.isTraceEnabled()) {
					log.trace("ADD NODE:\t" + newNode + " to " + currNode + "\t" + is_ins);
				}
			}
			if(currPtr.type() == PointerType.MIS) {
				this.addNodeAlignedToGraph(newNode, currNode);
				for(POGEdge inedge : new HashSet<POGEdge>(this.g.incomingEdgesOf(currNode))) {
					POGNode alignedPredecessor = this.g.getEdgeSource(inedge);
					if (!g.containsEdge(alignedPredecessor, newNode)) {
						POGEdge alternativeEdge = g.addEdge(alignedPredecessor, newNode);
						this.edgeReadIdx.put(alternativeEdge, new HashSet<Integer>());
					}
				}
			}
			prevNode = newNode;
		}

		if(log.isTraceEnabled()) {
			log.trace(
					(i + "\t" + j + "\t" 
					+ s1 + "\t" + s2 + "\t" 
					+ (is_ins ? "-" : nodeCH) + "\t"
					+ (is_del ? "-" : seqCH) + "\t" 
					+ currPos + "\t" + currPtr + "\t"
					+ (i != currPos.first()) + " "
					+ (j != currPos.second()) + "\t" 
					+ (i != currPos.first() && j != currPos.second()) + "\t" 
					+ prevNode+ "\t" 
					+ !(currPos.first() == 0 && currPos.second() == 0)));
		}
		
		if (currPos.equals(originPtr)) {
			return;
		}

		this.backtrackBuild(seq, 
							prevNode, 
							currPtr.ptr(),
							s1, 
							s2);
	}

	protected List<POGNode> backtrackPathRecurse(String seq,
												POGNode prevNode,
												Pair<Integer,Integer> currPos,
												String s1,
												String s2,
												List<POGNode> path) {

//		if (currPos.equals(originPtr)) {
		POGPtr currPtr = this.ptrMat.get(currPos.first(), currPos.second());
//		if (currPtr.ptr().equals(originPtr)) {
		if(currPos.equals(originPtr) || 
		   (currPtr.ptr().equals(originPtr) && currPtr.type() == PointerType.JUMP)) {
			this.queryIntval = new Interval(currPos.second(), -1);
			Collections.reverse(path);
			return path;
		}
		
//		POGPtr currPtr = this.ptrMat.get(currPos.first(), currPos.second());
		
		int i = currPos.first();
		int j = currPos.second();

		boolean is_del = currPtr.type() == PointerType.DEL;
		boolean is_ins = currPtr.type() == PointerType.INS;
		
		POGNode currNode = i > 0 ? nodeList.get(i-1) : null;
		char nodeCH = (i == 0 ? '-' : nodeList.get(i-1).getCH());
		char seqCH = (j == 0 ? '-' : seq.charAt(j-1));
		
		s1 += (is_ins ? "-" : nodeCH);
		s2 += (is_del ? "-" : seqCH);

		// if match, add index to node
		if(currPtr.type() == PointerType.MATCH) {
			prevNode = currNode;
			path.add(currNode);
		} else if(currPtr.type() == PointerType.DEL) {
			// don't add vertex, keep prevNode the same
		}
		else if(currPtr.type() == PointerType.INS) {
			path.add(null);
		}
		// if mismatch, add node and connect it to all prev nodes
		else {
			POGNode newNode = new POGNode(seqCH, this.readIDX);
			prevNode = newNode;
			path.add(null);
		}
		
//		System.out.println(
		if(log.isTraceEnabled())
			log.trace(
				(i + "\t" + j + "\t" 
				+ s1 + "\t" + s2 + "\t" 
				+ (is_ins ? "-" : nodeCH) + "\t"
				+ (is_del ? "-" : seqCH) + "\t" 
				+ currPos + "\t" + currPtr + "\t"
				+ (i != currPos.first()) + " "
				+ (j != currPos.second()) + "\t" 
				+ (i != currPos.first() && j != currPos.second()) + "\t" 
				+ prevNode+ "\t" 
				+ currPtr.ptr() + "\t"
				+ !(currPos.first() == 0 && currPos.second() == 0)));
		
		return this.backtrackPathRecurse(seq, 
										prevNode, 
										currPtr.ptr(),
										s1, 
										s2,
										path);
	}
	
	/**
	 * Backtrack and add nodes to POGraph
	 * @param seq
	 */
	protected void backtrack(String seq) {
		Pair<Integer,Integer> endCell = this.backtrackFromCell(seq.length());
		this.invalidateTopologyCache();
//		Pair<Integer,Integer> nextCell = this.alignPtrMat.get(endCell.first(), endCell.second());
//		POGPtr endPtr = this.backtrackFromPtr(seq.length());
//		System.out.println("END CELL:\t" + endCell);
//		POGPtr nextPtr = this.ptrMat.get(endPtr.ptr().first(), endPtr.ptr().second());
		
//		this.backtrackBuild(seq, null, nextPtr, endPtr, "", "");
//		this.backtrackBuild(seq, null, endPtr, null, "", "");
		this.backtrackBuild(seq, null, endCell, "", "");
	}
	
	/**
	 * Backtrack and build a path in the POGraph by recursing (tail-recursive).
	 * @param seq
	 * @param prevNode
	 * @param currP
	 * @param prevP
	 * @param s1
	 * @param s2
	 * @param path
	 * @return
	 * @deprecated
	 */
	protected List<POGNode> backtrackPathRecurse(String seq,
												 POGNode prevNode,
												 Pair<Integer,Integer> currP, 
												 Pair<Integer,Integer> prevP,
												 String s1,
												 String s2,
												 List<POGNode> path) {
		
		if(currP.equals(originPtr)) {
			this.queryIntval = new Interval(prevP.second(), -1);
			Collections.reverse(path);
			return path;
		}
		
		int i = prevP.first();
		int j = prevP.second();
		
		boolean is_del = prevP != null && (prevP.second() == currP.second()); 
		boolean is_ins = prevP != null && (prevP.first() == currP.first());
		
		POGNode currNode = i > 0 ? nodeList.get(i-1) : null;
		char nodeCH = (i == 0 ? '-' : nodeList.get(i-1).getCH());
		char seqCH = (j == 0 ? '-' : seq.charAt(j-1));
		
		s1 += (is_ins ? "-" : nodeCH);
		s2 += (is_del ? "-" : seqCH);
		
		// if match, add index to node
		if(nodeCH == seqCH) {
			prevNode = currNode;
			path.add(currNode);
		}
		else if(is_del) {
			// don't add vertex, keep prevNode the same
		}
		else if(is_ins) {
			path.add(null);
		}
		// if mismatch, add node and connect it to all prev nodes 
		else {				
			POGNode newNode = new POGNode(seqCH, this.readIDX);
			prevNode = newNode;
			path.add(null);
		}		

//		System.out.println(
		if(log.isTraceEnabled())
			log.trace((i + "\t" + j + "\t" 
					+ s1 + "\t" + s2 + "\t" 
					+ (is_ins ? "-" : nodeCH) + "\t" 
					+ (is_del ? "-" : seqCH) + "\t"
					+ currP + "\t" 
					+ prevP + "\t" 
					+ (i != currP.first()) + " " +  (j != currP.second()) + "\t" 
					+ (i != currP.first() && j != currP.second()) + "\t" 
					+ prevNode + "\t" 
					+ !(currP.first() == 0 && currP.second() == 0) + "\t"
					+ alignScoreMat.get(currP.first(), currP.second())));
		
		return this.backtrackPathRecurse(seq, prevNode, alignPtrMat.get(currP.first(), currP.second()), currP, s1, s2, path);
	}
	
	/**
	 * Backtracks a the sequence that has already been used to fill the matrix using {@link #alignSequence(String)}.
	 * @param seq
	 * @return
	 */
	protected List<POGNode> backtrackPath(String seq) {
		Pair<Integer,Integer> endCell = this.backtrackFromCell(seq.length());
		return this.backtrackPath(seq, endCell);
	}

	private List<POGNode> backtrackPath(String seq, Pair<Integer,Integer> endCell) {
		List<POGNode> path = this.backtrackPathRecurse(seq, null, endCell, "", "", new LinkedList<POGNode>());
		
		this.queryIntval = new Interval(this.queryIntval.start(), endCell.second());
		
		return path;
	}
	
	/**
	 * Add edges to alignedG (graph of aligned-to nodes)
	 */
	private void addNodeAlignedToGraph(POGNode newNode, POGNode algnNode) {
		if(!alignedG.containsVertex(algnNode)) return;
//		log.debug("ALIGNED TO:\t" + newNode + " => " + algnNode);
		alignedG.addVertex(newNode);
		alignedG.addEdge(algnNode, newNode);

		Set<POGEdge> incEdges = new HashSet<POGEdge>(alignedG.edgesOf(algnNode));
		for(POGEdge e : incEdges) {
			POGNode dest = alignedG.getEdgeTarget(e);
			if(dest == newNode) {
				dest = alignedG.getEdgeSource(e);
			}
			alignedG.addEdge(newNode, dest);
		}
	}
	
	/**
	 * Compile all nodes that are aligned to the given node
	 * @param {@code currNode} to get all other aligned to nodes
	 * @return {@code Set<POGNode>} all nodes aligned to the given node, including {@code currNode}
	 */
	private Set<POGNode> getAllAlignedToNodes(POGNode currNode) {
		Set<POGNode> alignedNodes = new HashSet<POGNode>();
		ArrayDeque<POGNode> pending = new ArrayDeque<POGNode>();
		alignedNodes.add(currNode);
		pending.add(currNode);
		while (!pending.isEmpty()) {
			POGNode node = pending.removeFirst();
			if (!alignedG.containsVertex(node)) {
				continue;
			}
			for (POGEdge edge : alignedG.edgesOf(node)) {
				POGNode source = alignedG.getEdgeSource(edge);
				POGNode target = alignedG.getEdgeTarget(edge);
				POGNode neighbor = node.equals(source) ? target : source;
				if (alignedNodes.add(neighbor)) {
					pending.addLast(neighbor);
				}
			}
		}
		return alignedNodes;
	}

	/** Build aligned components and topologically order the resulting quotient DAG. */
	private List<List<POGNode>> computeAlignmentColumns() {
		computeInternalDS();
		final List<List<POGNode>> components = new ArrayList<List<POGNode>>();
		final Map<POGNode,Integer> componentByNode = new HashMap<POGNode,Integer>();
		final List<Integer> componentOrderKeys = new ArrayList<Integer>();

		for (POGNode node : nodeList) {
			if (componentByNode.containsKey(node)) {
				continue;
			}
			Set<POGNode> connected = getAllAlignedToNodes(node);
			List<POGNode> component = new ArrayList<POGNode>(connected.size());
			int componentIndex = components.size();
			int orderKey = Integer.MAX_VALUE;
			for (POGNode alignedNode : connected) {
				if (!g.containsVertex(alignedNode)) {
					throw new IllegalStateException("aligned node is absent from the POA graph");
				}
				Integer previous = componentByNode.put(alignedNode, componentIndex);
				if (previous != null && previous.intValue() != componentIndex) {
					throw new IllegalStateException("node belongs to multiple alignment columns");
				}
				component.add(alignedNode);
				orderKey = Math.min(orderKey, alignedNode.getCounter());
			}
			components.add(component);
			componentOrderKeys.add(orderKey);
		}

		List<Set<Integer>> outgoing = new ArrayList<Set<Integer>>(components.size());
		int[] incomingCounts = new int[components.size()];
		for (int i = 0; i < components.size(); i++) {
			outgoing.add(new HashSet<Integer>());
		}
		for (POGEdge edge : g.edgeSet()) {
			int sourceColumn = componentByNode.get(g.getEdgeSource(edge));
			int targetColumn = componentByNode.get(g.getEdgeTarget(edge));
			if (sourceColumn != targetColumn
					&& outgoing.get(sourceColumn).add(targetColumn)) {
				incomingCounts[targetColumn]++;
			}
		}

		PriorityQueue<Integer> ready = new PriorityQueue<Integer>(Math.max(1, components.size()),
				new Comparator<Integer>() {
					@Override
					public int compare(Integer left, Integer right) {
						int keyComparison = Integer.compare(
								componentOrderKeys.get(left), componentOrderKeys.get(right));
						return keyComparison != 0 ? keyComparison : Integer.compare(left, right);
					}
				});
		for (int i = 0; i < incomingCounts.length; i++) {
			if (incomingCounts[i] == 0) {
				ready.add(i);
			}
		}

		List<List<POGNode>> ordered = new ArrayList<List<POGNode>>(components.size());
		while (!ready.isEmpty()) {
			int column = ready.remove();
			ordered.add(components.get(column));
			for (int target : outgoing.get(column)) {
				if (--incomingCounts[target] == 0) {
					ready.add(target);
				}
			}
		}
		if (ordered.size() != components.size()) {
			throw new IllegalStateException(
					"aligned-node columns form a cycle in the POA quotient graph");
		}
		return ordered;
	}
	
	/**
	 * Align a sequence to the graph
	 * @param seq {@code String} to align to the partial order graph
	 * @deprecated
	 */
	protected void alignSequenceOld(String seq) {
		computeInternalDS();
		
		TopologicalOrderIterator<POGNode,POGEdge> titer =
				new TopologicalOrderIterator<POGNode,POGEdge>(g);
		int i = 1;
		while(titer.hasNext()) {
			POGNode currN = titer.next();
			// get nodes of incoming edges 
			Set<POGEdge> adjE = g.incomingEdgesOf(currN);
			int numE = adjE.size();
			ArrayList<Integer> prevIDX = new ArrayList<Integer>(numE);
			ArrayList<POGNode> adjL = new ArrayList<POGNode>(numE);
			for(POGEdge e : adjE) {
				POGNode prevNode = g.getEdgeSource(e);
				adjL.add( prevNode );
				prevIDX.add( nodeMap.get(prevNode)+1 ); // force 1-index nodeIDX
			}
			
			for(int j = 1; j < seq.length()+1; j++) {
			
				char prevCh = seq.charAt(j-1);
//				System.out.print(
				if(log.isTraceEnabled())
					log.trace(
						(i) + "\t" 
							+ j + "\t"
							+ currN.toString() + "\t" 
							+ prevCh +"\t" 
							+ prevIDX + "\t" 
							+ adjL);
				fillCell(i, j, prevCh, currN.getCH(), prevIDX, adjL);
			}
			
			i++;
		}
		
//		log.debug("ADDED SEQ OF LEN " + seq.length() + "\t" 
//				+ this.g.vertexSet().size() + "\t"
//				+ this.g.edgeSet().size());		
	}
	
	protected void alignSequence(String seq) {
		computeInternalDS();
		prepareMatrices(seq.length());
		
		for(int i = 1; i <= this.nodeList.size(); i++) {
			int[] predecessors = this.predecessorRows[i - 1];
			char graphBase = this.nodeBases[i - 1];
			for(int j = 1; j < seq.length()+1; j++) {
				fillCell(i, j, seq.charAt(j-1), graphBase, predecessors);
			}
		}
	}

	/** Allocate matrices for exactly the graph and query currently being aligned. */
	protected void prepareMatrices(int sequenceLength) {
		int rows = this.nodeList.size() + 1;
		int columns = sequenceLength + 1;
		if (this.alignScoreMat instanceof ReusableAlignmentMatrixScore) {
			((ReusableAlignmentMatrixScore) this.alignScoreMat).resize(rows, columns);
		} else {
			this.alignScoreMat = new ReusableAlignmentMatrixScore(rows, columns);
		}
		this.coherenceScoreMat.resize(rows, columns);
		this.alignPtrMat.resize(rows, columns);
		this.ptrMat.resize(rows, columns);
		this.alignPtrMat.clearCell(0, 0);
		this.ptrMat.clearCell(0, 0);
		for (int row = 0; row < rows; row++) {
			this.coherenceScoreMat.set(row, 0, 0);
		}
		for (int column = 0; column < columns; column++) {
			this.coherenceScoreMat.set(0, column, 0);
		}
		this.initMatrix();
		this.afterMatrixInitialization();
	}

	/** Hook for alignment modes with graph-dependent boundary initialization. */
	protected void afterMatrixInitialization() {
		// Default boundary initialization is fully handled by initMatrix().
	}
	
	/**
	 * Fills a cell in the score/ptr matrix
	 * @param i
	 * @param j
	 * @param sequenceBase
	 * @param graphBase
	 */
	protected void fillCell(int i, int j, char sequenceBase, char graphBase,
			int[] predecessorRows) {
		int bestScore = 0;
		int bestCoherence = 0;
		int bestPreviousRow = 0;
		int bestPreviousColumn = 0;
		PointerType bestType = PointerType.UNDEF;
		boolean hasBest = false;
		int substitutionScore = scoreMat.score(graphBase, sequenceBase);
		PointerType substitutionType = graphBase == sequenceBase
				? PointerType.MATCH : PointerType.MIS;

		// Candidate order and strict-greater replacement preserve the old tie break.
		if (predecessorRows.length == 0) {
			bestScore = alignScoreMat.get(0, j - 1) + substitutionScore;
			bestCoherence = coherenceScoreMat.get(0, j - 1);
			bestPreviousColumn = j - 1;
			bestType = substitutionType;
			hasBest = true;
		}
		for (int predecessorIndex = 0;
				predecessorIndex < predecessorRows.length; predecessorIndex++) {
			int predecessorRow = predecessorRows[predecessorIndex];
			int score = alignScoreMat.get(predecessorRow, j - 1) + substitutionScore;
			int coherence = colorCoherenceTieBreaking
					? coherenceScoreMat.get(predecessorRow, j - 1)
							+ predecessorColorCoherence(i, predecessorRow) : 0;
			if (!hasBest || score > bestScore
					|| (score == bestScore && coherence > bestCoherence)) {
				bestScore = score;
				bestCoherence = coherence;
				bestPreviousRow = predecessorRow;
				bestPreviousColumn = j - 1;
				bestType = substitutionType;
				hasBest = true;
			}
		}
		for (int predecessorIndex = 0;
				predecessorIndex < predecessorRows.length; predecessorIndex++) {
			int predecessorRow = predecessorRows[predecessorIndex];
			int score = alignScoreMat.get(predecessorRow, j) + scoreMat.gapPenalty();
			int coherence = colorCoherenceTieBreaking
					? coherenceScoreMat.get(predecessorRow, j)
							+ predecessorColorCoherence(i, predecessorRow) : 0;
			if (score > bestScore || (score == bestScore && coherence > bestCoherence)) {
				bestScore = score;
				bestCoherence = coherence;
				bestPreviousRow = predecessorRow;
				bestPreviousColumn = j;
				bestType = PointerType.DEL;
			}
		}
		int insertionScore = alignScoreMat.get(i, j - 1) + scoreMat.gapPenalty();
		int insertionCoherence = colorCoherenceTieBreaking
				? coherenceScoreMat.get(i, j - 1) : 0;
		if (insertionScore > bestScore
				|| (insertionScore == bestScore && insertionCoherence > bestCoherence)) {
			bestScore = insertionScore;
			bestCoherence = insertionCoherence;
			bestPreviousRow = i;
			bestPreviousColumn = j - 1;
			bestType = PointerType.INS;
		}

		alignScoreMat.set(i, j, bestScore);
		coherenceScoreMat.set(i, j, bestCoherence);
		alignPtrMat.set(i, j, bestPreviousRow, bestPreviousColumn);
		ptrMat.set(i, j, bestPreviousRow, bestPreviousColumn, bestScore, bestType);
	}

	private int predecessorColorCoherence(int currentRow, int predecessorRow) {
		if (currentRow <= 0 || predecessorRow <= 0) {
			return 0;
		}
		int overlap = predecessorColorOverlap(currentRow, predecessorRow);
		int currentColorCount = nodeList.get(currentRow - 1).getNumReads();
		return currentColorCount == 0 ? 0
				: overlap * COLOR_COHERENCE_SCALE / currentColorCount;
	}

	private int predecessorColorOverlap(int currentRow, int predecessorRow) {
		int[] rows = predecessorRows[currentRow - 1];
		for (int index = 0; index < rows.length; index++) {
			if (rows[index] == predecessorRow) {
				return predecessorColorOverlaps[currentRow - 1][index];
			}
		}
		Set<Integer> currentColors =
				new HashSet<Integer>(nodeList.get(currentRow - 1).getIDX());
		int overlap = 0;
		for (int color : nodeList.get(predecessorRow - 1).getIDX()) {
			if (currentColors.contains(color)) {
				overlap++;
			}
		}
		return overlap;
	}

	/** Compatibility bridge for subclasses using the former boxed predecessor API. */
	@Deprecated
	protected void fillCell(int i, int j, char sequenceBase, char graphBase,
			List<Integer> predecessorIndices, List<POGNode> predecessorNodes) {
		int[] primitiveRows = new int[predecessorIndices.size()];
		for (int index = 0; index < primitiveRows.length; index++) {
			primitiveRows[index] = predecessorIndices.get(index);
		}
		fillCell(i, j, sequenceBase, graphBase, primitiveRows);
	}
	
	/**
	 * Compute the nodeMap (POGNode =&gt; index) and nodeList 
	 * (index =&gt; POGNode) data-structures
	 */
	protected void computeInternalDS() {
		if (topologyCacheValid
				&& cachedVertexCount == g.vertexSet().size()
				&& cachedEdgeCount == g.edgeSet().size()) {
			return;
		}
		nodeMap = this.computeOrder();
//		log.debug("NODE MAP:\t" + nodeMap);
		nodeList = new ArrayList<POGNode>(nodeMap.size());
		for(int i = 0; i < nodeMap.size(); i++) { nodeList.add(i, null); }		
		
		for(POGNode currNode : nodeMap.keySet()) {
			int index = nodeMap.get(currNode);
//			log.debug(index + "\t" + currNode);
			nodeList.set(index, currNode);
		}
		nodeBases = new char[nodeList.size()];
		predecessorRows = new int[nodeList.size()][];
		predecessorColorOverlaps = new int[nodeList.size()][];
		for (int index = 0; index < nodeList.size(); index++) {
			POGNode node = nodeList.get(index);
			nodeBases[index] = node.getCH();
			Set<POGEdge> incoming = g.incomingEdgesOf(node);
			int[] rows = new int[incoming.size()];
			int[] overlaps = new int[incoming.size()];
			Set<Integer> currentColors = new HashSet<Integer>(node.getIDX());
			int predecessorIndex = 0;
			for (POGEdge edge : incoming) {
				POGNode predecessor = g.getEdgeSource(edge);
				rows[predecessorIndex] = nodeMap.get(predecessor) + 1;
				int overlap = 0;
				for (int color : predecessor.getIDX()) {
					if (currentColors.contains(color)) {
						overlap++;
					}
				}
				overlaps[predecessorIndex] = overlap;
				predecessorIndex++;
			}
			predecessorRows[index] = rows;
			predecessorColorOverlaps[index] = overlaps;
		}
		cachedVertexCount = g.vertexSet().size();
		cachedEdgeCount = g.edgeSet().size();
		topologyCacheValid = true;
		topologyBuildCount++;
		
//		log.debug("NODE LIST:\t" + nodeList);
	}

	/** Exact row and predecessor order used by the dense DP implementation. */
	POAAlignmentTopologyView alignmentTopologyView() {
		computeInternalDS();
		List<POGNode> nodes = new ArrayList<POGNode>(nodeList);
		int[][] predecessors = new int[predecessorRows.length][];
		for (int row = 0; row < predecessorRows.length; row++) {
			predecessors[row] = predecessorRows[row].clone();
		}
		return new POAAlignmentTopologyView(nodes, predecessors);
	}

	protected void invalidateTopologyCache() {
		this.topologyCacheValid = false;
	}

	int topologyBuildCount() {
		return topologyBuildCount;
	}
	
	/**
	 * Compute the ordering of nodes, performed after the graph is perturbed.
	 * Runs topological sort and records the ordering.
	 */
	protected HashMap<POGNode,Integer> computeOrder() {
		HashMap<POGNode,Integer> nodeIDX = 
				new HashMap<POGNode,Integer>(g.vertexSet().size());
		TopologicalOrderIterator<POGNode,POGEdge> titer =
				new TopologicalOrderIterator<POGNode,POGEdge>(g);
		int i = 0;
		while(titer.hasNext()) {
			POGNode currN = titer.next();
			nodeIDX.put(currN, i);
			i++;
		}
		return nodeIDX;
	}
	
	/**
	 * 
	 * @return
	 */
	protected int getNumCols() {
		return computeAlignmentColumns().size();
	}	
	
	public int numNodes() {
		return this.g.vertexSet().size();
	}
	
	/**
	 * The {@link Interval} of the query sequence that was aligned last using
	 * {@link #alignToGraph(String)}
	 * @return
	 */
	public Interval queryInterval() {
		return this.queryIntval;
	}	
	
	/**
	 * 
	 * @return
	 */
	@Override
	public Iterator<POGNode> iterator() {
		return new POGIterator(this.g, this.alignedG);
	}
	
	private class POGIterator implements Iterator<POGNode> {
		private final Iterator<POGNode> columns;

		POGIterator(DefaultDirectedGraph<POGNode,POGEdge> pog,
				UndirectedGraph<POGNode,POGEdge> ag) {
			List<POGNode> representatives = new ArrayList<POGNode>();
			for (List<POGNode> column : AbstractPartialOrderGraph.this.computeAlignmentColumns()) {
				representatives.add(column.get(0));
			}
			this.columns = representatives.iterator();
		}

		public boolean hasNext() {
			return columns.hasNext();
		}

		public POGNode next() {
			if (!hasNext()) {
				throw new NoSuchElementException();
			}
			return columns.next();
		}
	}
	
	
}
