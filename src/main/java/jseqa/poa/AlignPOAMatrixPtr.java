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
	private int[] colPtrMat;
	
	/**
	 * an integer matrix for the row pointer
	 */
	private int[] rowPtrMat;

	private int capacityRows;
	private int capacityColumns;
	
	public AlignPOAMatrixPtr(int nrows, int ncols) {
		super(1, 1);
		resize(nrows, ncols);
	}

	public void resize(int nrows, int ncols) {
		if (rowPtrMat == null || nrows > capacityRows || ncols > capacityColumns) {
			capacityRows = grow(capacityRows, nrows);
			capacityColumns = grow(capacityColumns, ncols);
			rowPtrMat = new int[capacityRows * capacityColumns];
			colPtrMat = new int[capacityRows * capacityColumns];
		}
		this.numRows = nrows;
		this.numCols = ncols;
	}

	@Override
	public void clearCell(int i, int j) {
		colPtrMat[index(i, j)] = -1;
		rowPtrMat[index(i, j)] = -1;
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
		rowPtrMat[index(i, j)] = iPtr;
		colPtrMat[index(i, j)] = jPtr;
	}
	
	/**
	 * Get the pointer pair from cell [i][j]
	 * @param i
	 * @param j
	 * @return
	 */
	public Pair<Integer,Integer> get(int i, int j) {
		return new Pair<Integer,Integer>(row(i, j), column(i, j));
	}

	public int row(int i, int j) {
		return rowPtrMat[index(i, j)];
	}

	public int column(int i, int j) {
		return colPtrMat[index(i, j)];
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
		return "[" + row(i, j) + ";" + column(i, j) + "]";
	}

	private int index(int row, int column) {
		return row * capacityColumns + column;
	}

	private static int grow(int current, int required) {
		int grown = current == 0 ? 1 : current;
		while (grown < required) {
			grown = grown < 1024 ? grown * 2 : grown + grown / 2;
		}
		return grown;
	}
}
