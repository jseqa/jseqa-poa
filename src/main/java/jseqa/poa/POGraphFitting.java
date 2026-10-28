package jseqa.poa;

import java.util.Collection;
import jseqa.alignment.substitution.SubstitutionMatrix;
import jseqa.poa.POGPtr.PointerType;
import jseqa.util.Pair;

/**
 * Partial order graph that aligns a sequence to the graph using fitting alignment, 
 * fitting the entire input sequence to the existing graph. This likely behaves the most 
 * as you would expect.
 *
 */
public class POGraphFitting extends AbstractPartialOrderGraph {

	public POGraphFitting() {
		super();
	}

	public POGraphFitting(POGNodeFactory nodeFactory) {
		super(nodeFactory);
	}
	
	/**
	 * {@inheritDoc}
	 */
	public POGraphFitting(String initseq) {
		super(initseq);
	}

	public POGraphFitting(String initseq, POGNodeFactory nodeFactory) {
		super(initseq, nodeFactory);
	}

	public POGraphFitting(Collection<String> alignment) {
		super(alignment);
	}

	public POGraphFitting(Collection<String> alignment, POGNodeFactory nodeFactory) {
		super(alignment, nodeFactory);
	}
	
	/**
	 * {@inheritDoc}
	 */
	public POGraphFitting(SubstitutionMatrix sm) {
		super(sm);
	}

	public POGraphFitting(SubstitutionMatrix sm, POGNodeFactory nodeFactory) {
		super(sm, nodeFactory);
	}
	
	/**
	 * {@inheritDoc}
	 */
	public POGraphFitting(String initSeq, SubstitutionMatrix sm) {
		super(initSeq, sm);
	}

	public POGraphFitting(String initSeq, SubstitutionMatrix sm, POGNodeFactory nodeFactory) {
		super(initSeq, sm, nodeFactory);
	}

	public POGraphFitting(Collection<String> alignment, SubstitutionMatrix sm) {
		super(alignment, sm);
	}

	public POGraphFitting(Collection<String> alignment, SubstitutionMatrix sm,
			POGNodeFactory nodeFactory) {
		super(alignment, sm, nodeFactory);
	}
	
	/**
	 * {@inheritDoc}
	 */
	public POGraphFitting(AbstractPartialOrderGraph pog) {
		super(pog);
	}
	
	/**
	 * {@inheritDoc}
	 * <br>
	 * Provide free-rides to the start of the sequence, i.e., start at any node in the POGraph
	 */
	@Override
	protected void initMatrix() {
		// penalize not using the sequence
		for(int j = 0; j < alignScoreMat.ncols(); j++) {
			alignScoreMat.set(0, j, j*scoreMat.gapPenalty());
			if(j >= 1) { 
				alignPtrMat.set(0, j, 0, j-1);
				this.ptrMat.set(0, j, new POGPtr(new Pair<Integer,Integer>(0,j-1), alignScoreMat.get(0, j), PointerType.INS));
			}
		}
		// reward using the entire sequence
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
	 * <br>
	 * Selects backtracking cell amongst the best scoring from the maximum column (i.e., the full sequence).
	 */
	@Override
	protected Pair<Integer, Integer> backtrackFromCell(int maxSeqLen) {
		Pair<Integer,Integer> maxCell = new Pair<Integer,Integer>(-1,-1);
		double maxScore = -10000000.0;
		
		int maxCol = maxSeqLen;
		for(int i = 0; i < this.numNodes()+1; i++) {
			if(this.alignScoreMat.get(i, maxCol) > maxScore) {
				maxScore = this.alignScoreMat.get(i, maxCol);
				maxCell = new Pair<Integer,Integer> (i, maxCol);
			}
		}
		return maxCell;
	}

	@Override
	@Deprecated
	protected POGPtr backtrackFromPtr(int maxSeqLen) {
		// TODO Auto-generated method stub
		return null;
	}

}
