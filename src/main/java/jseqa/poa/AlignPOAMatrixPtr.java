package jseqa.poa;

import jseqa.alignment.AbstractAlignmentMatrix;
import jseqa.util.Pair;

/**
 * Alignment pointer matrix for the {@link com.digitalproteomics.sequence.poa.POGraph}
 *
 */
public class AlignPOAMatrixPtr extends AbstractAlignmentMatrix {

	/**
	 * an integer matrix for the column pointer
	 */
	private int[][] colPtrMat;
	
	/**
	 * an integer matrix for the row pointer
	 */
	private int[][] rowPtrMat;
	
	public AlignPOAMatrixPtr(int nrows, int ncols) {
		super(nrows, ncols);
		colPtrMat = new int[nrows][ncols];
		rowPtrMat = new int[nrows][ncols];
	}

	@Override
	public void clearCell(int i, int j) {
		colPtrMat[i][j] = -1;
		rowPtrMat[i][j] = -1;
	}

	@Override
//	public void bottomOutCell(int i, int j) {
	public void setLowestVal(int i, int j) {
		clearCell(i, j);
	}
	
	/**
	 * Set the pointer of cell [i,j] to [iPtr, jPtr]
	 * @param i
	 * @param j
	 * @param iPtr
	 * @param jPtr
	 */
	public void set(int i, int j, int iPtr, int jPtr) {
		rowPtrMat[i][j] = iPtr;
		colPtrMat[i][j] = jPtr;
	}
	
	/**
	 * Get the pointer pair from cell [i][j]
	 * @param i
	 * @param j
	 * @return
	 */
	public Pair<Integer,Integer> get(int i, int j) {
		return new Pair<Integer,Integer>(rowPtrMat[i][j], colPtrMat[i][j]);
	}
	
//	/**
//	 * overload to output matrix as a String
//	 */
//	@Override 
//	public String toString() {
//		String str = "";
//		for(int i = 0; i < this.numRows; i++) {
//			for(int j = 0; j < this.numCols; j++) {
//				str += "[" + rowPtrMat[i][j] + ";" + colPtrMat[i][j] + "]";
//				if(j < this.numCols-1) str += ",";
//			}
//			str += "\n";
//		}
//		return str;
//	}

	@Override
	public String cellToString(int i, int j) {
		return "[" + rowPtrMat[i][j] + ";" + colPtrMat[i][j] + "]";
	}	
}
