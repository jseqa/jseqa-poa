package jseqa.poa;

import java.util.Collection;
import java.util.List;
import java.util.Set;

import jseqa.util.Interval;
import org.jgrapht.graph.DefaultDirectedGraph;

/**
 * Interface for all partial order alignment graphs.
 *
 */
public interface PartialOrderAlignment {

	/**
	 * Adds a sequence to the partial order alignment graph
	 * @param seq {@code String} sequence to add to the graph
	 */
	public void add(String seq);
	
	/**
	 * Adds all sequences in the collection to the partial order alignment graph.
	 * @param seqs {@code Collection} of {@code String} sequences to add
	 */
	public void addAll(Collection<String> seqs);
	
	/**
	 * Align the sequence to the POGraph without adding nodes/edges, and return the backtracked list of 
	 * {@link com.digitalproteomics.sequence.poa.POGNode} representing the path of the aligned sequence.
	 * @param seq {@code String} sequence to align to the POGraph
	 * @return List of {@code POGNodes} representing the path of the aligned sequence
	 */
	public List<POGNode> alignToGraph(String seq);
	
	/**
	 * Compute the multiple sequence alignment (MSA) from the partial order graph.
	 * @return {@code List<String>} of each row of the MSA
	 */
	public List<String> getAlignment();
	
	/**
	 * Compute the consensus sequence, as represented by the heaviest path
	 * @return {@code String} consensus sequence of the heaviest path.
	 */
	public String consensus();
	
	/**
	 * Output graph to a DOT file.
	 * @param dotFilename {@code String} DOT filename
	 */
	public void toDOT(String dotFilename);

	public void toDOT(String dotFilename, boolean showNodeLabels);

	/**
	 * Get the internal graph representation of the POGraph for analysis.
	 * @return
	 */
	public DefaultDirectedGraph<POGNode,POGEdge> getGraph();

	/**
	 * Get the number of nodes in the partial order alignment graph
	 * @return
	 */
	public int numNodes();
	
	/**
	 * The {@link jseqa.util.digitalproteomics.sequence.poa.Interval} of the query sequence that was aligned last using
	 * {@link #alignToGraph(String)}
	 * @return
	 */
	public Interval queryInterval();

	public Set<Integer> getEdgeIdx(POGEdge e);
}
