package jseqa.poa;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.List;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.config.Configuration;
import org.apache.logging.log4j.core.config.LoggerConfig;
import org.junit.Test;

import jseqa.alignment.substitution.DnaSubstitutionMatrix;

public class POGraphGlobalTest {
	
	@Test
	public void testMore() {
		List<String> seqs = Arrays.asList(
			"AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA",
			"TGTACCTGTTTGTGAGGCTA",
			"AGTTCCTGCTGCGTTTGCTGGACTTATGTACTTGTTTGTGAGGCAA",
		   "AAGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGGTTTGTGGAGGCAA",
			"AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTTAGGCAA",
			"AGTTCCTGCTGCGTTTGCT",
			"AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTT",			
			"AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA"
			);
		//  AAGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGGTTTGTGGAGGCAA
		String exp = "AAGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGGTTTGTGGAGGCAA";
		POGraphGlobal pog = new POGraphGlobal(seqs.get(0), new DnaSubstitutionMatrix(-15));
//		POGraphGlobal pog = new POGraphGlobal(seqs.get(0));
		
//		System.out.println(pog.startNodes.iterator().next().getFirst().getCounter());
//		System.out.println(pog.endNodes.iterator().next().getFirst().getCounter());

		LoggerContext ctx = (LoggerContext) LogManager.getContext(false);
		Configuration config = ctx.getConfiguration();
		LoggerConfig loggerConfig = config.getLoggerConfig(LogManager.ROOT_LOGGER_NAME); 
//		loggerConfig.setLevel(Level.TRACE);
		loggerConfig.setLevel(Level.DEBUG);
		ctx.updateLoggers();
		
		for(int i = 1; i < seqs.size(); i++) {
//		for(int i = 1; i < 3; i++) {
			pog.add(seqs.get(i));
			System.out.println(pog.startNodes);
			System.out.println(pog.endNodes);
			System.out.println("-------------------------");
		}
		pog.toDOT("/tmp/test_graph2.dot");
		
		String cons = pog.consensus();
		List<String> algn = pog.getAlignment();
		for(String s : algn) {
			System.out.println(s);
		}
		
		for(int i = 0; i < seqs.size(); i++) {
//		for(int i = 6; i < 7; i++) {
			List<POGNode> path = pog.alignToGraph(seqs.get(i));
			for(int ni = 0; ni < pog.numNodes(); ni++) {
//				for(int nj = 0; nj < seqs.get(i).length()+1; nj++) {
				for(int nj = 34; nj < seqs.get(i).length()+1; nj++) {
//					System.err.println(ni + "\t" + nj + "\t" + pog.alignScoreMat.get(ni, nj) + "\t" + pog.ptrMat.get(ni, nj));
				}
			}
			
//			List<POGNode> readPath = pog.getPathOfReadIdx(i);
//			System.out.println("READ PATH:\t" + readPath);
//			for(int pos = 0; pos < readPath.size(); pos++) {
//				System.out.println(i + "\t" + pos + "\t" + seqs.get(i).charAt(pos) + "\t" + readPath.get(pos).getCH() + "\t" + readPath.get(pos));
//				assertEquals(seqs.get(i).charAt(pos), readPath.get(pos).getCH());
//			}
			
			int pos = 0;
			for(POGNode cnode : path) {
//				System.out.println(i + "\t" + cnode + "\t" + pos + "\t" + seqs.get(i).charAt(pos));
//				assertTrue(cnode.getIDX().contains(i));
//				assertEquals(seqs.get(i).charAt(pos), cnode.getCH());
				pos++;
			}
		}
	}

}
