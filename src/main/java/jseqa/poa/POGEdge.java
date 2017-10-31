package jseqa.poa;

import org.jgrapht.graph.DefaultEdge;

/**
 * Edge for a partial order alignment graph
 *
 */
public class POGEdge extends DefaultEdge {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	private POGEdgeType et;
	
	/**
	 * 
	 */
	public POGEdge() {
		et = POGEdgeType.ALIGN_TO;
	}
	
	/**
	 * 
	 * @param e
	 */
	public POGEdge(POGEdgeType e) {
		et = e;
	}
	
	/**
	 * 
	 * @return
	 */
	public POGEdgeType getEdgeType() {
		return et;
	}	
}
