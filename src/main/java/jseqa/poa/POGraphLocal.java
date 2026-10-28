package jseqa.poa;

import java.util.Collection;
import java.util.List;

import jseqa.alignment.substitution.SubstitutionMatrix;
import jseqa.poa.POGPtr.PointerType;
import jseqa.util.Pair;

/**
 * Partial order graph that performs local alignment of sequences to the graph.
 *
 */
public class POGraphLocal extends AbstractPartialOrderGraph {
	
	public POGraphLocal() {
		super();
	}

	public POGraphLocal(POGNodeFactory nodeFactory) {
		super(nodeFactory);
	}
		
	/**
	 * {@inheritDoc}
	 */
	public POGraphLocal(String initseq) {
		super(initseq);
	}

	public POGraphLocal(String initseq, POGNodeFactory nodeFactory) {
		super(initseq, nodeFactory);
	}

	public POGraphLocal(Collection<String> alignment) {
		super(alignment);
	}

	public POGraphLocal(Collection<String> alignment, POGNodeFactory nodeFactory) {
		super(alignment, nodeFactory);
	}
	
	/**
	 * {@inheritDoc}
	 */
	public POGraphLocal(SubstitutionMatrix sm) {
		super(sm);
	}

	public POGraphLocal(SubstitutionMatrix sm, POGNodeFactory nodeFactory) {
		super(sm, nodeFactory);
	}
	
	/**
	 * {@inheritDoc}
	 */
	public POGraphLocal(String initSeq, SubstitutionMatrix sm) {
		super(initSeq, sm);
	}

	public POGraphLocal(String initSeq, SubstitutionMatrix sm, POGNodeFactory nodeFactory) {
		super(initSeq, sm, nodeFactory);
	}

	public POGraphLocal(Collection<String> alignment, SubstitutionMatrix sm) {
		super(alignment, sm);
	}

	public POGraphLocal(Collection<String> alignment, SubstitutionMatrix sm,
			POGNodeFactory nodeFactory) {
		super(alignment, sm, nodeFactory);
	}
	
	/**
	 * {@inheritDoc}
	 */
	public POGraphLocal(AbstractPartialOrderGraph pog) {
		super(pog);
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void fillCell(int i, 
			   int j, 
			   char prevSeqCh,
			   char prevGraphCh,
			   List<Integer> predIDX, 
			   List<POGNode> predN) {
		super.fillCell(i, j, prevSeqCh, prevGraphCh, predIDX, predN);
				
		if(this.alignScoreMat.get(i, j) <= 0) {
			this.alignScoreMat.set(i, j, 0);
			this.alignPtrMat.set(i, j, 0, 0);
			this.ptrMat.set(i, j, new POGPtr(new Pair<Integer,Integer>(0,0), 0, PointerType.JUMP));
		}
	}
	
	/**
	 * Initializes the score matrix and backtrack pointer matrix, defaulting
	 * to free-rides on the first column/row
	 */
	@Override
	protected void initMatrix() {
		for(int j = 0; j < alignScoreMat.ncols(); j++) {
			alignScoreMat.set(0, j, 0);			
			if(j >= 1) { 
				alignPtrMat.set(0, j, 0, 0);
				this.ptrMat.set(0, j, new POGPtr(new Pair<Integer,Integer>(0,0), 0, PointerType.JUMP));
			}
		}
		for(int i = 0; i < alignScoreMat.nrows(); i++) {
			alignScoreMat.set(i, 0, 0);
			if( i >= 1) {
				alignPtrMat.set(i, 0, 0, 0);
				this.ptrMat.set(i, 0, new POGPtr(new Pair<Integer,Integer>(0,0), 0, PointerType.JUMP));
			}
		}
	}
	
	/**
	 * {@inheritDoc}
	 * 
	 * Backtracks from maximum scoring cell from entire matrix, as per local alignment.
	 */
	@Override
	protected Pair<Integer,Integer> backtrackFromCell(int maxSeqLen) {
		Pair<Integer,Integer> maxCell = new Pair<Integer,Integer>(-1,-1);
		double maxScore = -10000000.0;
		for(int i = 0; i <= this.numNodes(); i++) {
			for(int j = 0; j <= maxSeqLen; j++) {
				if(this.alignScoreMat.get(i, j) > maxScore) {
					maxScore = this.alignScoreMat.get(i, j);
					maxCell = new Pair<Integer,Integer> (i, j);
				}
			}
		}		
		return maxCell;
	}

	@Deprecated
	@Override
	protected POGPtr backtrackFromPtr(int maxSeqLen) {
		Pair<Integer,Integer> endcell = this.backtrackFromCell(maxSeqLen);
		return this.ptrMat.get(endcell.first(), endcell.second());
	}
}
