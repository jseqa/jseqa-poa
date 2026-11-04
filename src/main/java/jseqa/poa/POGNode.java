package jseqa.poa;

import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Node for a partial order alignment graph
 * Contains:
 * <ul>
 * <li>Char (nt/aa)</li>
 * <li>Read indices</li>
 * </ul>
 * Each object has a unique identifier, so no two objects will hash the same.
 */
public class POGNode {
	
	/**
	 * char representing node
	 */
	private char ch;
	
	/**
	 * List of index (i.e., indexing the sequences)
	 */
	private List<Integer> idx;
	
	/**
	 * counter of number of instances
	 */
	private static final AtomicInteger INSTANCE_COUNTER = new AtomicInteger(-1);
	
	/**
	 * instance counter accessible to the object
	 */
	private int counter;
	
	/**
	 * Constructor that creates an empty list and undefined char (N)
	 */
	public POGNode() {
		ch = 'N';
		init();
	}
	
	/**
	 * Constructor that assigns the char and index
	 * @param n char to assign to the node
	 * @param nidx int index to add
	 */
	public POGNode(char n, int nidx) {
		ch = n;
		init();
		idx.add(nidx);
	}

	/** Create an independent node with the same residue and read membership. */
	POGNode(POGNode other) {
		ch = other.ch;
		init();
		idx.addAll(other.idx);
	}
	
	/**
	 * 
	 */
	private void init() {
		idx = new LinkedList<Integer>();
		counter = INSTANCE_COUNTER.incrementAndGet();
	}
	
	/**
	 * add a new read index to list
	 * @param newIDX
	 */
	public void addIDX(int newIDX) {
		idx.add(newIDX);
	}
	
	/**
	 * get the {@code List<Integer>} of indices at the node
	 * @return
	 */
	public List<Integer> getIDX() {
		return idx;
	}
	
	/**
	 * Determines if node is a dummy node, i.e., if only has an ID of -1
	 * @return
	 */
	public boolean isDummyNode() {
		return idx.size() == 1 && idx.contains(-1);
	}
	
	/**
	 * Get the number of reads stored the node, if is dummy node, returns 0
	 * @return int number of read indices stored in node.
	 */
	public int getNumReads() {
		return this.isDummyNode() ? 0 : idx.size();
	}
	
	/**
	 * get the character at the node
	 * @return char 
	 */
	public char getCH() {
		return ch;
	}
	
	/**
	 * 
	 */
	@Override
	public String toString() {
		return "[" + this.counter + ": " + this.ch + " " + this.idx.toString() + "]";
	}
	
	/**
	 * Get the counter of the object, the object ID.
	 * @return int object index identifier
	 */
	public int getCounter() {
		return counter;
	}
	
	/**
	 * 
	 */
	@Override
	public int hashCode() {
		int val = 7;
		val += 31 + counter;
		return val;
	}
	
	/**
	 * 
	 */
	public boolean equals(Object obj) {
		if(obj == null) { return false; }
		if(obj instanceof POGNode) {
			POGNode pog = (POGNode)obj;
			return (pog.counter == counter);
		}
		return false;
	}
}
