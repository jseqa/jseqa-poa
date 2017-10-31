package jseqa.poa;

import jseqa.util.Pair;

/**
 * Represents a pointer in the partial order graph, used in computing the 
 * maximum value in the recurrence relation.
 *
 */
public class POGPtr implements Comparable<POGPtr> {

	private Pair<Integer,Integer> ptr;
	
	private int score;
	
	private PointerType ptrType;
	
	/**
	 * Default constructor
	 */
	public POGPtr() {
		this(new Pair<Integer,Integer>(-1,-1), 0, PointerType.UNDEF);
	}
	
	/**
	 * Constructor that defines the location of the pointer, score, and type
	 * @param p
	 * @param val int value
	 * @param pt {@code PointerType} type of pointer
	 */
	public POGPtr(Pair<Integer,Integer> p, int val, PointerType pt) {
		score  = val;
		ptr = p;
		ptrType = pt;
	}
	
	public POGPtr(Pair<Integer,Integer> p, int val) {
		this(p, val, PointerType.UNDEF);
	}
	
	@Deprecated
	public int getScore() {
		return score;
	}
	
	@Deprecated
	public Pair<Integer,Integer> getPtr() {
		return ptr;
	}
	
	@Deprecated
	public PointerType getType() {
		return ptrType;
	}
	
	public int score() {
		return score;
	}
	
	public Pair<Integer,Integer> ptr() {
		return ptr;
	}
	
	public PointerType type() {
		return ptrType;
	}
	
	@Override
	public int compareTo(POGPtr o) {
		return Integer.compare(score, o.getScore());
	}
	
	public String toString() {
		return "[" + ptr.toString() + "; " + score + "; " + ptrType + "]";
	}
	
	/**
	 * Represents the different types of pointers
	 *
	 */
	public enum PointerType {
		/** Match */
		MATCH,
		/** Insertion */
		INS,
		/** Deletion */
		DEL,
		/** Mismatch */
		MIS,
		/** Jump */
		JUMP,
		/** Undefined */
		UNDEF;
		
		public String toString() {
			if(this == PointerType.MATCH) { return "M"; }
			else if(this == PointerType.INS) { return "I"; }
			else if(this == PointerType.DEL) { return "D"; }
			else if(this == PointerType.MIS) { return "X"; }
			else if(this == PointerType.JUMP) { return "J"; }
			else { return "U"; }
		}
	}
}
