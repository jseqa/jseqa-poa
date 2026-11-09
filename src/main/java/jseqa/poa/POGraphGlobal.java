package jseqa.poa;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import jseqa.alignment.substitution.SubstitutionMatrix;
import jseqa.poa.POGPtr.PointerType;
import jseqa.util.Pair;

/**
 * Partial order graph that performs global alignment of sequences to the graph. Matches 
 * to best start/end node in the graph.
 *
 */
public class POGraphGlobal extends AbstractPartialOrderGraph {

	/** Set of start nodes */
	protected Map<POGNode,Integer> startNodes;
	
	/** Set of end nodes */
	protected Map<POGNode,Integer> endNodes;
	
	public POGraphGlobal() {
		super();
		this.initNodeSets();
	}

	public POGraphGlobal(POGNodeFactory nodeFactory) {
		super(nodeFactory);
		this.initNodeSets();
	}
		
	/** {@inheritDoc} */
	public POGraphGlobal(String initseq) {
		super(initseq);
		this.initNodeSets();
		this.initStartNodesInMatrix();
	}

	public POGraphGlobal(String initseq, POGNodeFactory nodeFactory) {
		super(initseq, nodeFactory);
		this.initNodeSets();
		this.initStartNodesInMatrix();
	}

	public POGraphGlobal(Collection<String> alignment) {
		super(alignment);
		this.initNodeSets();
		this.initializeEndpointMapsFromGraph();
	}

	public POGraphGlobal(Collection<String> alignment, POGNodeFactory nodeFactory) {
		super(alignment, nodeFactory);
		this.initNodeSets();
		this.initializeEndpointMapsFromGraph();
	}
	
	/** {@inheritDoc} */
	public POGraphGlobal(SubstitutionMatrix sm) {
		super(sm);
		this.initNodeSets();
	}

	public POGraphGlobal(SubstitutionMatrix sm, POGNodeFactory nodeFactory) {
		super(sm, nodeFactory);
		this.initNodeSets();
	}
	
	/** {@inheritDoc} */
	public POGraphGlobal(String initSeq, SubstitutionMatrix sm) {
		super(initSeq, sm);
		this.initNodeSets();
		this.initStartNodesInMatrix();
	}

	public POGraphGlobal(String initSeq, SubstitutionMatrix sm, POGNodeFactory nodeFactory) {
		super(initSeq, sm, nodeFactory);
		this.initNodeSets();
		this.initStartNodesInMatrix();
	}

	public POGraphGlobal(Collection<String> alignment, SubstitutionMatrix sm) {
		super(alignment, sm);
		this.initNodeSets();
		this.initializeEndpointMapsFromGraph();
	}

	public POGraphGlobal(Collection<String> alignment, SubstitutionMatrix sm,
			POGNodeFactory nodeFactory) {
		super(alignment, sm, nodeFactory);
		this.initNodeSets();
		this.initializeEndpointMapsFromGraph();
	}
	
	/** {@inheritDoc} */
	public POGraphGlobal(AbstractPartialOrderGraph pog) {
		super(pog);
		this.initNodeSets();
		this.initializeEndpointMapsFromGraph();
	}
	
	private void initNodeSets() {
		this.startNodes = new HashMap<POGNode,Integer>();
		this.endNodes = new HashMap<POGNode,Integer>();
	}
	
	/** {@inheritDoc} */
	@Override
	public void add(String seq) {
		super.add(seq);
		this.computeInternalDS();
		this.addStartEndNodes();
	}
	
	@Override
	protected void initGraph(String seq) {
		super.initGraph(seq);
		super.computeInternalDS();
	}
	
	/**
	 * Go through all nodes, add first/last of current read to start/end set
	 */
	private void addStartEndNodes() {
		int currReadIdx = this.readIDX-1;
		POGNode firstCurrIdxNode = null;
		POGNode lastCurrIdxNode = null;
		for(int i = 0; i < this.nodeList.size(); i++) {
			POGNode currNode = this.nodeList.get(i);
			if(currNode.getIDX().contains(currReadIdx)) {
				if(firstCurrIdxNode == null) {
					firstCurrIdxNode = currNode;
				}
				lastCurrIdxNode = currNode;
			}
		}
		if (firstCurrIdxNode != null) {
			this.startNodes.put(firstCurrIdxNode, this.nodeList.indexOf(firstCurrIdxNode) + 1);
			this.endNodes.put(lastCurrIdxNode, this.nodeList.indexOf(lastCurrIdxNode) + 1);
		}
	}
	
	/**
	 * Updates the index of the start/end nodes
	 */
	protected void updateStartEndNodeIndex() {
		for(int i = 0; i < this.nodeList.size(); i++) {
			POGNode currNode = this.nodeList.get(i);
			if(this.startNodes.containsKey(currNode)) {
				this.startNodes.put(currNode, i + 1);
			}
			if(this.endNodes.containsKey(currNode)) {
				this.endNodes.put(currNode, i + 1);
			}
		}
	}
	
	/**
	 * Determine existing start/end nodes in the POGraph and initialize all start nodes in the matrix of the 
	 * POGraph initial row to have no penalty.
	 */
	protected void initStartNodesInMatrix() {
		if(this.g.vertexSet().size() == 0) {
			return;
		}
		if (this.startNodes.isEmpty() || this.endNodes.isEmpty()) {
			this.initializeEndpointMapsFromGraph();
		}
		
		for(int i = 0; i < alignScoreMat.nrows(); i++) {
			alignScoreMat.set(i, 0, i*scoreMat.gapPenalty());
			if( i >= 1 && i <= this.nodeList.size()) {
				POGNode currNode = this.nodeList.get(i-1);
				if(this.startNodes.containsKey(currNode)) {
					this.alignScoreMat.set(i, 0, 0);
					alignPtrMat.set(i, 0, 0, 0);
					this.ptrMat.set(i, 0, 0, 0, 0, PointerType.JUMP);
				}
			}
		}
	}

	private void initializeEndpointMapsFromGraph() {
		for (int i = 0; i < this.nodeList.size(); i++) {
			POGNode node = this.nodeList.get(i);
			if (g.inDegreeOf(node) == 0) {
				this.startNodes.put(node, i + 1);
			}
			if (g.outDegreeOf(node) == 0) {
				this.endNodes.put(node, i + 1);
			}
		}
	}

	@Override
	protected void afterMatrixInitialization() {
		this.initStartNodesInMatrix();
	}
	
	/**
	 * {@inheritDoc}
	 * <br>
	 * Penalizes linearly the first row/column, but does not penalize any of the start 
	 * nodes (performed in {@linkplain #initStartNodesInMatrix()}.
	 */
	protected void initMatrix() {
		for(int j = 0; j < alignScoreMat.ncols(); j++) {
			alignScoreMat.set(0, j, j*scoreMat.gapPenalty());
			if(j >= 1) { 
				alignPtrMat.set(0, j, 0, j-1);
				this.ptrMat.set(0, j, 0, j-1, alignScoreMat.get(0, j), PointerType.INS);
			}
		}
		for(int i = 0; i < alignScoreMat.nrows(); i++) {
			alignScoreMat.set(i, 0, i*scoreMat.gapPenalty());
			if( i >= 1) {
				alignPtrMat.set(i, 0, i-1, 0);
				this.ptrMat.set(i, 0, i-1, 0, alignScoreMat.get(i, 0), PointerType.DEL);
			}
		}
	}

	/**
	 * Return the index of the node to start backtracking from, index from the score matrix (not the node list index).
	 * @param max_seq_len int max sequence index
	 * @return int index of node to start backtracking from
	 */
	protected int initBacktrackFromNodeIndex(int max_seq_len) {
		int maxScore = -1;
		int maxIDX = -1;
		// nodeList is 0-index, but score matrix is 1-index
		for(int i = 1; i < this.nodeList.size()+1; i++) {
			POGNode currNode = this.nodeList.get(i-1);
			int oDeg = this.g.outDegreeOf(currNode);

			int currScore = this.alignScoreMat.get(i, max_seq_len);
			if(log.isTraceEnabled() && oDeg == 0) {
				log.trace("END:\t" + i + "\t" + oDeg + "\t" 
							+ currScore + " > " + maxScore + "\t" 
							+ this.nodeList.get(i-1) + "\t" 
							+ this.ptrMat.get(i-1, max_seq_len));
			}
			if(oDeg == 0 && currScore > maxScore) {
				maxIDX = i; //i+1;
				maxScore = currScore;
			}
		}
		return maxIDX >= 0 ? maxIDX : this.nodeList.size();
	}
	
	/**
	 * {@inheritDoc}
	 * 
	 * select the best score among all {@linkplain #endNodes}.
	 */
	protected Pair<Integer,Integer> backtrackFromCell(int maxSeqLen) {
//		return new Pair<Integer,Integer>(this.initBacktrackFromNodeIndex(maxSeqLen), maxSeqLen);
		this.updateStartEndNodeIndex();
		Pair<Integer,Integer> maxCell = new Pair<Integer,Integer>(-1,-1);
		for(POGNode endNode : this.endNodes.keySet()) {
			int currEndIdx = this.endNodes.get(endNode);
			if(this.isBetterAlignmentCell(currEndIdx, maxSeqLen,
					maxCell.first(), maxCell.second())) {
				maxCell = new Pair<Integer,Integer>(currEndIdx, maxSeqLen);
			}
		}
		return maxCell;
	}

	/** Endpoint row order used by dense single-best tie handling. */
	int[] alignmentEndRowsInSelectionOrder() {
		this.updateStartEndNodeIndex();
		int[] rows = new int[this.endNodes.size()];
		int index = 0;
		for (POGNode node : this.endNodes.keySet()) {
			rows[index++] = this.endNodes.get(node);
		}
		return rows;
	}

	@Override
	protected boolean isKBestEndCell(int row, int column, int queryLength) {
		if (column != queryLength) {
			return false;
		}
		this.updateStartEndNodeIndex();
		return this.endNodes.containsValue(row);
	}
	
	@Override
	protected POGPtr backtrackFromPtr(int maxSeqLen) {
		Pair<Integer,Integer> endcell = this.backtrackFromCell(maxSeqLen);
		return new POGPtr(endcell, 0, PointerType.UNDEF);
	}
}
