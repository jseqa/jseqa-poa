package jseqa.poa;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.List;

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
		POGraphGlobal pog = new POGraphGlobal(seqs.get(0), new DnaSubstitutionMatrix(-15));
//		POGraphGlobal pog = new POGraphGlobal(seqs.get(0));
		
//		System.out.println(pog.startNodes.iterator().next().getFirst().getCounter());
//		System.out.println(pog.endNodes.iterator().next().getFirst().getCounter());

		for(int i = 1; i < seqs.size(); i++) {
//		for(int i = 1; i < 3; i++) {
			pog.add(seqs.get(i));
		}
		String cons = pog.consensus();
		List<String> algn = pog.getAlignment();
		assertEquals(seqs.size(), algn.size());
		assertFalse(cons.isEmpty());
		
		for(int i = 0; i < seqs.size(); i++) {
//		for(int i = 6; i < 7; i++) {
			List<POGNode> path = pog.alignToGraph(seqs.get(i));
			assertFalse(path.isEmpty());
//			List<POGNode> readPath = pog.getPathOfReadIdx(i);
//			System.out.println("READ PATH:\t" + readPath);
//			for(int pos = 0; pos < readPath.size(); pos++) {
//				System.out.println(i + "\t" + pos + "\t" + seqs.get(i).charAt(pos) + "\t" + readPath.get(pos).getCH() + "\t" + readPath.get(pos));
//				assertEquals(seqs.get(i).charAt(pos), readPath.get(pos).getCH());
//			}
			
		}
	}

}
