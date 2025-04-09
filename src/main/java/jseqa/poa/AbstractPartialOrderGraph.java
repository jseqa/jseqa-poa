package jseqa.poa;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Optional;
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
 * Provides the skeletal implementation of Partial Order Graphs
 *
 */
public abstract class AbstractPartialOrderGraph implements Iterable<POGNode>, PartialOrderAlignment {

	final static Logger log = LogManager.getLogger(AbstractPartialOrderGraph.class);
	
	final private int MAX_NODES = 2000;
	
	final private int MAX_SEQ_LEN = 500;
	
	final private Pair<Integer,Integer> originPtr = new Pair<Integer,Integer>(0,0);
	
	/** Index of sequences added to POGraph */
	protected int readIDX;
	
	/** Map from POGNode =&gt; index */
	private HashMap<POGNode,Integer> nodeMap;
	
	/** List of nodes, providing a map back from index =&gt; POGNode */
	protected ArrayList<POGNode> nodeList;
	
	/** scoring matrix, rows =&gt; graph, cols =&gt; sequence */
//	protected AlignMatrixScore alignScoreMat;
	protected AlignmentMatrixScore alignScoreMat;
	
	/**
	 * matrix of pointers to cells
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
	
	/** {@link #Interval} of the aligned query sequence */
	protected Interval queryIntval;

	Map<POGEdge, Set<Integer>> edgeReadIdx;

	/**
	 * Default constructor.
	 */
	public AbstractPartialOrderGraph() {
		this.init();
		this.initMatrix();
	}
	
	/**
	 * Constructor that initializes with a single sequence.
	 * @param initseq {@code String} sequence to initialize.
	 */
	public AbstractPartialOrderGraph(String initseq) {
		this.init();
		this.initMatrix();
		this.initGraph(initseq);
	}
	
	/**
	 * Constructor that specifies a particular scoring matrix to use
	 * @param sm {@link #ScoreMatrix} to use for scoring the alignments
	 */
	public AbstractPartialOrderGraph(SubstitutionMatrix sm) {
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
		this.init(sm);
		this.initMatrix();
		this.initGraph(initSeq);
	}
	
	/**
	 * Copy constructor that performs a shallow copy of graph data-strucutres
	 * @param otherPOG
	 */
	public AbstractPartialOrderGraph(AbstractPartialOrderGraph otherPOG) {
		alignScoreMat = new AlignmentMatrixScore(MAX_NODES, MAX_SEQ_LEN);
		alignPtrMat = new AlignPOAMatrixPtr(MAX_NODES, MAX_SEQ_LEN);
		this.ptrMat = new AlignPOGPtrMatrix(MAX_NODES, MAX_SEQ_LEN);
		
		this.g = (DefaultDirectedGraph<POGNode, POGEdge>) otherPOG.g.clone();
		this.alignedG = otherPOG.alignedG;
		
		scoreMat = otherPOG.scoreMat;
		readIDX = otherPOG.readIDX;
		this.edgeReadIdx = new HashMap<>(otherPOG.edgeReadIdx);
		this.computeInternalDS();
		this.initMatrix();
	}
	
	/**
	 * Adds a sequence to the partial order alignment graph
	 * @param seq {@code String} sequence to add to the graph
	 */
	public void add(String seq) {
		if(this.g.vertexSet().size() == 0) {
			this.initGraph(seq);
		} else {
			this.alignSequence(seq);
			this.backtrack(seq);
			readIDX++;
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
	 * Align the sequence to the POGraph without adding nodes/edges, and return the backtracked list of 
	 * {@link #POGNodes} representing the path of the aligned sequence.
	 * @param seq {@code String} sequence to align to the POGraph
	 * @return List of {@code POGNodes} representing the path of the aligned sequence
	 */
	public List<POGNode> alignToGraph(String seq) {
		this.alignSequence(seq);
		return this.backtrackPath(seq);
	}
	
	/**
	 * Compute the multiple sequence alignment (MSA) from the partial order graph.
	 * @return {@code List<String>} of each row of the MSA
	 */
	public List<String> getAlignment() {
		computeInternalDS();
		
		int numNodes = g.vertexSet().size();
		int numCols = this.getNumCols()+10;
//		log.debug("NODES:\t" + numNodes + "\tCOLS:\t" + numCols);
		
		int numReads = this.readIDX;
		List<StringBuilder> algn = new ArrayList<StringBuilder>(numReads);
		for(int i = 0; i < numReads; i++) {
			StringBuilder sb = new StringBuilder(numCols);
			for(int j = 0; j < numCols; j++) {
				sb.append('-');
			}
			algn.add(sb);
		}
		
		Set<POGNode> visitedNodes = new HashSet<POGNode>(alignedG.vertexSet().size());
		int index = 0;
		for(int i = 0; i < numNodes; i++) {
			POGNode currNode = nodeList.get(i);
			
//			log.debug(i + "\t" 
//					+ index + "\t" 
//					+ currNode + "\t" 
//					+ visitedNodes.contains(currNode));
			
			if(visitedNodes.contains(currNode)) { continue; }
			Set<POGNode> alignedToNodes = this.getAllAlignedToNodes(currNode);
			
			for(POGNode alignedNode : alignedToNodes) {
				char nodeCH = alignedNode.getCH();
				List<Integer> nodeIDX = alignedNode.getIDX();
				for(Integer idx : nodeIDX) {
					if(idx < 0) { continue; }
					algn.get(idx).setCharAt(index, nodeCH);
				}
			}
			// add to visited nodes
			visitedNodes.addAll(alignedToNodes);
			index++;
		}
		
		ArrayList<String> alignLst = new ArrayList<String>(numReads);
		for(int i = 0; i < numReads; i++) {
			alignLst.add(i, algn.get(i).toString());
		}
		
		return alignLst;
	}

	public DefaultDirectedGraph<POGNode,POGEdge> getGraph() {
		return this.g;
	}

	/**
	 * Compute the consensus sequence, as represented by the heaviest path
	 * @return {@code String} consensus sequence of the heaviest path.
	 */
	public String consensus() {
		
		TopologicalOrderIterator<POGNode,POGEdge> titer = 
				new TopologicalOrderIterator<POGNode,POGEdge>(g);
		HashMap<POGNode,Integer> nodeScores = new HashMap<POGNode,Integer>(g.vertexSet().size());
		HashMap<POGNode,POGNode> nodePtrs = new HashMap<POGNode,POGNode>(g.vertexSet().size());
//		POGNode initNode = titer.next();
		POGNode initNode = null;
		nodeScores.put(initNode, -1);
		nodePtrs.put(initNode, null);
		int i = 1;		
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
				nodeScores.put(currN, 0);
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
		String consensusSeq = "";
		POGNode prevNode = maxEndNode;
		while(prevNode != null) {
			consensusSeq = prevNode.getCH() + consensusSeq;
//			log.debug("BACKTRACK:\t" + prevNode + "\t" + consensusSeq);
			prevNode = nodePtrs.get(prevNode);
		}
		
		return consensusSeq;
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
	 * {@link #DNAScoreMatrix} by as the default
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
		alignScoreMat = new AlignmentMatrixScore(MAX_NODES, MAX_SEQ_LEN);
		alignPtrMat = new AlignPOAMatrixPtr(MAX_NODES, MAX_SEQ_LEN);
		this.ptrMat = new AlignPOGPtrMatrix(MAX_NODES, MAX_SEQ_LEN);
		scoreMat = sm;
		readIDX = 0;
		this.edgeReadIdx = new HashMap<>();
	}
	
	/**
	 * initialize the graph with a single sequence, making a chain
	 * @param seq {@code String} initial sequence
	 */
	protected void initGraph(String seq) {
		POGNode prevNode = null;
		for(int i = 0; i < seq.length(); i++) {
			POGNode currNode =  new POGNode(seq.charAt(i), 0);
			g.addVertex( currNode );
			if(prevNode != null) {
				POGEdge e = g.addEdge(prevNode, currNode);
				this.safeEdgeAddIdx(e, readIDX);
			}
			alignedG.addVertex(currNode);
			prevNode = currNode;
		}	
		readIDX++;
	}	

	protected void safeEdgeAddIdx(POGEdge e, int idx) {
		if (!this.edgeReadIdx.containsKey(e)) {
			this.edgeReadIdx.put(e, new HashSet<>());
		}
		this.edgeReadIdx.get(e).add(idx);
	}

	/**
	 * Return the cell in the score matrix to start backtracking from, with first => node IDX, and second => seq IDX.
	 * @param maxSeqLen 
	 * @return
	 */
	abstract protected Pair<Integer,Integer> backtrackFromCell(int maxSeqLen);
	
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
			currNode.addIDX(readIDX);
//			log.trace("MAKE EDGE?:\t" + (prevNode != null) +  !g.containsEdge(currNode, prevNode));
			if (prevNode != null && !g.containsEdge(currNode, prevNode)) {
				g.addEdge(currNode, prevNode);
			}
			POGEdge e = g.getEdge(currNode, prevNode);
			this.safeEdgeAddIdx(e, readIDX);

			prevNode = currNode;
		} else if(currPtr.type() == PointerType.DEL) {
			// don't add vertex, keep prevNode the same
		} else if(currPtr.type() == PointerType.INS) {
			POGNode newNode = new POGNode(seqCH, this.readIDX);
			g.addVertex(newNode);
			if(prevNode != null) {
				g.addEdge(newNode, prevNode);
			}
			POGEdge e = g.getEdge(newNode, prevNode);
			this.safeEdgeAddIdx(e, readIDX);
			if(log.isTraceEnabled()) {
				log.trace("INSERTION:\t" + newNode + " <- " + prevNode);
			}
			prevNode = newNode;
		}
		// if mismatch, add node and connect it to all prev nodes
		else {
			POGNode newNode = new POGNode(seqCH, this.readIDX);
			g.addVertex(newNode);
			if(prevNode != null) {
				POGEdge e = g.addEdge(newNode, prevNode);
				this.safeEdgeAddIdx(e, readIDX);
				// add node/edges for alignedTo graph
//				if (!is_ins) {
				if(currPtr.type() == PointerType.MIS) { // prevP
					this.addNodeAlignedToGraph(newNode, currNode);
					// connect to all incoming edges
					for(POGEdge inedge : this.g.incomingEdgesOf(currNode)) {
						POGNode algnToPrevNode = this.g.getEdgeSource(inedge);
						POGEdge lclE = g.addEdge(algnToPrevNode, newNode);
						this.safeEdgeAddIdx(lclE, readIDX);
					}
				}
				if(log.isTraceEnabled()) {
					log.trace("ADD NODE:\t" + newNode + " to " + currNode + "\t" + is_ins);
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

		Set<POGEdge> incEdges = alignedG.edgesOf(algnNode);
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
		if(!alignedG.containsVertex(currNode)) { 
			Set<POGNode> singleSet = new HashSet<POGNode>();
			singleSet.add(currNode);
			return singleSet;
		}
		Set<POGEdge> incEdges = alignedG.edgesOf(currNode);
		Set<POGNode> alignedNodes = new HashSet<POGNode>(incEdges.size()+1);		
		alignedNodes.add(currNode);
		for(POGEdge e : incEdges) {
			POGNode src = alignedG.getEdgeSource(e);
			POGNode dest = alignedG.getEdgeTarget(e);
			alignedNodes.add(src);
			alignedNodes.add(dest);
		}
		return alignedNodes;
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
		
		for(int i = 1; i < this.nodeList.size()+1; i++) {
			POGNode currN = this.nodeList.get(i-1);
			
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
//				System.out.println(
				if(log.isTraceEnabled())
					log.trace(
							i + "\t" 
							+ j + "\t"
							+ currN.toString() + "\t" 
							+ prevCh +"\t" 
							+ prevIDX + "\t" 
							+ adjL);
				fillCell(i, j, prevCh, currN.getCH(), prevIDX, adjL);
			}
		}
	}
	
	/**
	 * Fills a cell in the score/ptr matrix
	 * @param i
	 * @param j
	 * @param prevSeqCh
	 * @param prevGraphCh
	 * @param predIDX
	 * @param predN
	 */
	protected void fillCell(int i, 
							int j, 
							char prevSeqCh, 
							char prevGraphCh,
							List<Integer> predIDX, 
							List<POGNode> predN) {
		ArrayList<POGPtr> scoreV = new ArrayList<POGPtr>(predIDX.size()*2+1);
		// starts
		if(predIDX.isEmpty()) {
			int val = this.scoreMat.score(prevGraphCh, prevSeqCh);
			scoreV.add(new POGPtr(new Pair<Integer,Integer>(i-1, j-1), val, PointerType.MATCH));
		}
		// diagonals
		for(int ii = 0; ii < predIDX.size(); ii++) {
			char a = predN.get(ii).getCH();
			int idx = predIDX.get(ii);
			int val = alignScoreMat.get(idx, j-1) 
					+ scoreMat.score(prevGraphCh, prevSeqCh);
			PointerType ptype = (prevGraphCh == prevSeqCh ? PointerType.MATCH : PointerType.MIS);
			scoreV.add(new POGPtr(new Pair<Integer,Integer>(idx, j-1), val, ptype)) ;
//			System.out.print(
//			log.trace("\t{"+ idx + " " 
//					+ val + "=" + alignScoreMat.get(idx, j-1) 
//					+ "+" + scoreMat.getScore(prevGraphCh, prevSeqCh));
		}
		// deletions
		for(int ii = 0; ii < predIDX.size(); ii++) {
			char a = predN.get(ii).getCH();
			int idx = predIDX.get(ii);
			int val = alignScoreMat.get(idx, j) 
					+ scoreMat.gapPenalty();
			scoreV.add(new POGPtr(new Pair<Integer,Integer>(idx, j), val, PointerType.DEL)) ;
		}
		// insertion
		int insScore = alignScoreMat.get(i, j-1)
						+ scoreMat.gapPenalty();
		scoreV.add(new POGPtr(new Pair<Integer,Integer>(i, j-1), insScore, PointerType.INS));
		
		// find max
		Optional<POGPtr> maxv = scoreV.stream()
									  .max(Comparator.comparing(p -> p.getScore()));
		POGPtr maxPOG = maxv.get();
		
		int maxScore = maxPOG.score();
		
//		log.trace("\t" + maxScore +"\t" +maxPOG + "\t" + scoreV);
		
		alignScoreMat.set(i, j, maxScore);
		alignPtrMat.set(i,  j, maxPOG.ptr().first(), maxPOG.ptr().second());
		this.ptrMat.set(i, j, maxPOG);
	}
	
	/**
	 * Compute the nodeMap (POGNode =&gt; index) and nodeList 
	 * (index =&gt; POGNode) data-structures
	 */
	protected void computeInternalDS() {
		nodeMap = this.computeOrder();
//		log.debug("NODE MAP:\t" + nodeMap);
		nodeList = new ArrayList<POGNode>(nodeMap.size());
		for(int i = 0; i < nodeMap.size(); i++) { nodeList.add(i, null); }		
		
		for(POGNode currNode : nodeMap.keySet()) {
			int index = nodeMap.get(currNode);
//			log.debug(index + "\t" + currNode);
			nodeList.set(index, currNode);
		}
		
//		log.debug("NODE LIST:\t" + nodeList);
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
		int count = 0;
		for(POGNode currNode : this) {
			count++;
		}
		return count;
	}	
	
	public int numNodes() {
		return this.g.vertexSet().size();
	}
	
	/**
	 * The {@link #Interval} of the query sequence that was aligned last using
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
        
        private TopologicalOrderIterator<POGNode,POGEdge> currNodeIter;
        private HashSet<POGNode> visitedNodes;
        private POGNode currNode;
        
        public POGIterator(DefaultDirectedGraph<POGNode,POGEdge> pog,
        				   UndirectedGraph<POGNode,POGEdge> ag) {
            currNodeIter = new TopologicalOrderIterator<POGNode,POGEdge>(pog);
//            currNode = currNodeIter.next();
            currNode = null;
            visitedNodes = new HashSet<POGNode>(pog.vertexSet().size());
            Set<POGNode> adj = this.getAllAlignedToNodes(currNode);
            visitedNodes.addAll(adj);
        }
        
        private Set<POGNode> getAllAlignedToNodes(POGNode currNode) {
    		if(!alignedG.containsVertex(currNode)) { 
    			Set<POGNode> singleSet = new HashSet<POGNode>();
    			singleSet.add(currNode);
    			return singleSet;
    		}
    		Set<POGEdge> incEdges = alignedG.edgesOf(currNode);
    		Set<POGNode> alignedNodes = new HashSet<POGNode>(incEdges.size()+1);		
    		alignedNodes.add(currNode);
    		for(POGEdge e : incEdges) {
    			POGNode src = alignedG.getEdgeSource(e);
    			POGNode dest = alignedG.getEdgeTarget(e);
    			alignedNodes.add(src);
    			alignedNodes.add(dest);
    		}
    		return alignedNodes;
    	}
        
        public boolean hasNext() {
            return currNodeIter.hasNext();
        }

        public POGNode next() {
        	if(this.hasNext()) {
        		Set<POGNode> alignedToNodes = this.getAllAlignedToNodes(currNode);
        		visitedNodes.addAll(alignedToNodes);
        		while(visitedNodes.contains(currNode)) {
        			currNode = currNodeIter.next();
        		}
        		return currNode;
        	}
        	throw new NoSuchElementException();
        }        
    }
	
	
}
