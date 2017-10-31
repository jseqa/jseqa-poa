package jseqa.poa;

import jseqa.alignment.AbstractAlignmentMatrix;

/**
 * Matrix holding backtracking pointers for a {@link com.digitalproteomics.sequence.poa.POGraph}.
 * An N x M matrix of {@link com.digitalproteomics.sequence.poa.POGPtr} objects used to backtrack
 * alignments of the partial order graph.
 *
 */
public class AlignPOGPtrMatrix extends AbstractAlignmentMatrix {	
	
	/** An integer matrix for the row pointer */
	private POGPtr[][] ptrMat;
	
	/**
	 * Constructor that initializes an {@code nrows} x {@code ncols} matrix of 
	 * backtracking pointers
	 * @param nrows int number of rows
	 * @param ncols int number of columns
	 */
	public AlignPOGPtrMatrix(int nrows, int ncols) {
		super(nrows, ncols);
		ptrMat = new POGPtr[nrows][ncols];
		for(int i = 0; i < nrows; i++) {
			for(int j = 0; j < ncols; j++) {
				ptrMat[i][j] = new POGPtr();
			}
		}
	}
	
	/** 
	 * {@inheritDoc} 
	 */
	@Override
	public void clearCell(int i, int j) {		
		ptrMat[i][j] = null;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
//	public void bottomOutCell(int i, int j) {
	public void setLowestVal(int i, int j) {
		clearCell(i, j);
	}
	
	/**
	 * Set the pointer of cell [i,j] to [iPtr, jPtr]
	 * @param i int row index
	 * @param j int column index
	 * @param ptr {@code POGPtr} backtracking pointer
	 */
	public void set(int i, int j, POGPtr ptr) {
		ptrMat[i][j] = ptr;
	}
	
	/**
	 * Get the pointer pair from cell [i][j]
	 * @param i int row index
	 * @param j int column index
	 * @return {@code POGPtr} at the cell [i][j]
	 */
	public POGPtr get(int i, int j) {
		return ptrMat[i][j];
	}
	
//	/**
//	 * overload to output matrix as a String
//	 */
//	@Override public String toString() {
//		String str = "";
//		for(int i = 0; i < this.numRows; i++) {
//			for(int j = 0; j < this.numCols; j++) {
//				str += ptrMat[i][j] + " ";
//				if(j < this.numCols-1) str += ",";
//			}
//			str += "\n";
//		}
//		return str;
//	}

	@Override
	public String cellToString(int i, int j) {
		return "" + this.ptrMat[i][j];
	}

}
