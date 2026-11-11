package jseqa.poa;

import java.util.Arrays;
import jseqa.alignment.AbstractAlignmentMatrix;
import jseqa.poa.POGPtr.PointerType;
import jseqa.util.Pair;

/**
 * Matrix holding backtracking pointers for a POGraph.
 * An N x M matrix of {@link jseqa.poa.POGPtr} objects used to backtrack
 * alignments of the partial order graph.
 *
 */
public class AlignPOGPtrMatrix extends AbstractAlignmentMatrix {	
	private static final PointerType[] POINTER_TYPES = PointerType.values();
	
	/** An integer matrix for the row pointer */
	private int[] previousRows;
	private int[] previousColumns;
	private int[] scores;
	private byte[] operations;
	private int capacityRows;
	private int capacityColumns;
	
	/**
	 * Constructor that initializes an {@code nrows} x {@code ncols} matrix of 
	 * backtracking pointers
	 * @param nrows int number of rows
	 * @param ncols int number of columns
	 */
	public AlignPOGPtrMatrix(int nrows, int ncols) {
		super(1, 1);
		resize(nrows, ncols);
	}

	public void resize(int nrows, int ncols) {
		if (previousRows == null || nrows > capacityRows || ncols > capacityColumns) {
			capacityRows = grow(capacityRows, nrows);
			capacityColumns = grow(capacityColumns, ncols);
			int size = capacityRows * capacityColumns;
			previousRows = new int[size];
			previousColumns = new int[size];
			scores = new int[size];
			operations = new byte[size];
			Arrays.fill(previousRows, -1);
			Arrays.fill(previousColumns, -1);
			Arrays.fill(operations, (byte) PointerType.UNDEF.ordinal());
		}
		this.numRows = nrows;
		this.numCols = ncols;
	}
	
	/** 
	 * {@inheritDoc} 
	 */
	@Override
	public void clearCell(int i, int j) {		
		set(i, j, -1, -1, 0, PointerType.UNDEF);
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
		set(i, j, ptr.ptr().first(), ptr.ptr().second(), ptr.score(), ptr.type());
	}

	public void set(int i, int j, int previousRow, int previousColumn,
			int score, PointerType type) {
		int index = index(i, j);
		previousRows[index] = previousRow;
		previousColumns[index] = previousColumn;
		scores[index] = score;
		operations[index] = (byte) type.ordinal();
	}
	
	/**
	 * Get the pointer pair from cell [i][j]
	 * @param i int row index
	 * @param j int column index
	 * @return {@code POGPtr} at the cell [i][j]
	 */
	public POGPtr get(int i, int j) {
		return new POGPtr(new Pair<Integer,Integer>(previousRow(i, j), previousColumn(i, j)),
				score(i, j), type(i, j));
	}

	public int previousRow(int i, int j) {
		return previousRows[index(i, j)];
	}

	public int previousColumn(int i, int j) {
		return previousColumns[index(i, j)];
	}

	public int score(int i, int j) {
		return scores[index(i, j)];
	}

	public PointerType type(int i, int j) {
		return POINTER_TYPES[operations[index(i, j)]];
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
		return "" + get(i, j);
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
