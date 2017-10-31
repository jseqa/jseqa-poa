package jseqa.poa;
/**
 * Represents the different types of edges in a POG.
 *
 */
public enum POGEdgeType {
	/**
	 * denotes a predecessor/successor link
	 */
	PRED_SUCC,
	/**
	 * denotes a connection between two aligned nodes
	 */
	ALIGN_TO;
	
	@Override
	public String toString() {
		if(super.toString() == "PRED_SUCC") { return "L"; }
		else if(super.toString() == "ALIGN_TO") { return "A"; }
		else { return "?"; }
	}
}
