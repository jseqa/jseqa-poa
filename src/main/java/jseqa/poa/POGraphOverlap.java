package jseqa.poa;

import java.util.Collection;
import jseqa.alignment.substitution.SubstitutionMatrix;
import jseqa.poa.POGPtr.PointerType;
import jseqa.util.Pair;

/**
 * Partial order alignment graph that does not penalize starts/ends (i.e., overlap alignment)
 *
 */
public class POGraphOverlap extends AbstractPartialOrderGraph {

	/**
	 * Default constructor.
	 */
	public POGraphOverlap() {
		super();
	}

	public POGraphOverlap(POGNodeFactory nodeFactory) {
		super(nodeFactory);
	}
	
	/**
	 * {@inheritDoc}
	 */
	public POGraphOverlap(String initseq) {
		super(initseq);
	}

	public POGraphOverlap(String initseq, POGNodeFactory nodeFactory) {
		super(initseq, nodeFactory);
	}

	public POGraphOverlap(Collection<String> alignment) {
		super(alignment);
	}

	public POGraphOverlap(Collection<String> alignment, POGNodeFactory nodeFactory) {
		super(alignment, nodeFactory);
	}
	
	/**
	 * {@inheritDoc}
	 */
	public POGraphOverlap(SubstitutionMatrix sm) {
		super(sm);
	}

	public POGraphOverlap(SubstitutionMatrix sm, POGNodeFactory nodeFactory) {
		super(sm, nodeFactory);
	}
	
	/**
	 * {@inheritDoc}
	 */
	public POGraphOverlap(String initSeq, SubstitutionMatrix sm) {
		super(initSeq, sm);
	}

	public POGraphOverlap(String initSeq, SubstitutionMatrix sm, POGNodeFactory nodeFactory) {
		super(initSeq, sm, nodeFactory);
	}

	public POGraphOverlap(Collection<String> alignment, SubstitutionMatrix sm) {
		super(alignment, sm);
	}

	public POGraphOverlap(Collection<String> alignment, SubstitutionMatrix sm,
			POGNodeFactory nodeFactory) {
		super(alignment, sm, nodeFactory);
	}

	/**
	 * {@inheritDoc}
	 */
	public POGraphOverlap(AbstractPartialOrderGraph pog) {
		super(pog);
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void initMatrix() {
		for(int j = 0; j < this.alignScoreMat.ncols(); j++) {
			this.alignScoreMat.set(0, j, 0);
			if(j >= 1) { 
				this.alignPtrMat.set(0, j, 0, 0);
				this.ptrMat.set(0, j, new POGPtr(new Pair<Integer,Integer>(0,0), 0, PointerType.JUMP));
			}
		}
		for(int i = 0; i < this.alignScoreMat.nrows(); i++) {
			this.alignScoreMat.set(i, 0, 0);
			if( i >= 1) {
				this.alignPtrMat.set(i, 0, 0, 0);
				this.ptrMat.set(i, 0, new POGPtr(new Pair<Integer,Integer>(0,0), 0, PointerType.JUMP));
			}
		}
	}
	
	/**
	 * {@inheritDoc}
	 * <br>
	 * The {@code POGraphOverlap} version searches for best value in last row/column
	 */
	@Override
	protected Pair<Integer,Integer> backtrackFromCell(int maxSeqLen) {
		Pair<Integer,Integer> maxCell = new Pair<Integer,Integer>(-1,-1);
		double maxScore = -10000000.0;
		
		int maxCol = maxSeqLen;
		for(int i = 0; i < this.numNodes()+1; i++) {
			if(this.alignScoreMat.get(i, maxCol) > maxScore) {
				maxScore = this.alignScoreMat.get(i, maxCol);
				maxCell = new Pair<Integer,Integer> (i, maxCol);
			}
		}
		int maxRow = this.numNodes();
		for(int j = 0; j < maxSeqLen+1; j++) {
			if(this.alignScoreMat.get(maxRow, j) > maxScore) {
				maxScore = this.alignScoreMat.get(maxRow, j);
				maxCell = new Pair<Integer,Integer> (maxRow, j);
			}
		}
		return maxCell;
	}

	@Override
	protected POGPtr backtrackFromPtr(int maxSeqLen) {
		Pair<Integer,Integer> endcell = this.backtrackFromCell(maxSeqLen);
		return this.ptrMat.get(endcell.first(), endcell.second());
	}
	
}
