package jseqa.poa;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.config.Configuration;
import org.apache.logging.log4j.core.config.LoggerConfig;
import org.junit.Test;

import jseqa.alignment.substitution.CustomDnaSubstitutionMatrix;

public class POGraphFittingTest {
	
	@SuppressWarnings("unchecked")
	@Test
	public void testJExample() {
		List<String> refSeqs = new ArrayList<String>() {{
			add("CTACTGGTACCTGGATCCCTGGGGCACAGGCACCCTGGTCACCATCTCTTCAG"); 	//IGHJ1*01
			add("TGATGCTTTTGATCCCTGGGGCCCAGGCACCCTGGTCACCGTCTCCTCAG");    	//IGHJ2*01
			add("TGATGCTTTCGATCCCTGGGGCCCAGGCACCCTGGTCACCGTCTCCTCAG");		//IGHJ2*02
			add("TGACACTCGGTTGGATCTCTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG");	//IGHJ3*01
			add("ACGATCGGTTGGATCCCTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG");		//IGHJ3*02
			add("ACTACTTTAACTTGTGGGGCCCAGGCACCCTGGTCACCGTCTCCTCAG");		//IGHJ4*01
			add("ACTACTTTAACATCTGGGGCCCAGGCACCCTGGTCACCGTCTCCTCAG");		//IGHJ4*02
			add("ATGACTGGTTGGATCTCTGGGGCCAGGGCACCCTGGTCACCGTCTCTTCA"); //N");		//IGHJ5*01
			add("ATGACTGGTTGGATCTCTGGGGCCAGGGCACCCTGGTCACTGTCTCTTCAG");		//IGHJ5*02
			add("ATTACTACGGCATGGACCTCTGGGGCCCAGGGACCCTCGTCACCGTCTCTTCAG");	//IGHJ6*01
			add("ATTACTACGGCATGGACCCCTGGGGCCCAGGGACCCTCGTCACCGTCTCTTCAG");	//IGHJ6*02
		}};
		List<String> ids = new ArrayList<String>() {{
			add("IGHJ1*01");
			add("IGHJ2*01");
			add("IGHJ2*02");
			add("IGHJ3*01");
			add("IGHJ3*02");
			add("IGHJ4*01");
			add("IGHJ4*02");
			add("IGHJ5*01");
			add("IGHJ6*01");
			add("IGHJ6*02");
		}};
		
		CustomDnaSubstitutionMatrix scoremat = new CustomDnaSubstitutionMatrix(5, -4, -4, -1);

		POGraphFitting poaG = new POGraphFitting(scoremat);
		for(String s : refSeqs) {
			poaG.add(s);
		}
		
		for(int i = 0; i < refSeqs.size(); i++) {
			String seq = refSeqs.get(i);
			List<POGNode> path = poaG.alignToGraph(seq);
			for(POGNode currNode : path) {
				assertTrue(currNode.getIDX().contains(i));
			}
		}
	}
	
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
		POGraphFitting pog = new POGraphFitting(seqs.get(0));

		for(int i = 1; i < seqs.size(); i++) {
			pog.add(seqs.get(i));
		}
		List<String> aln = pog.getAlignment();
		
		int num = 0;
		for(POGNode currNode : pog) {
			num++;
		}
		assertEquals(50, num);
		
		for(int i = 0; i < seqs.size(); i++) {
			List<POGNode> path = pog.alignToGraph(seqs.get(i));
			int pos = 0;
			for(POGNode cnode : path) {
				assertTrue(cnode.getIDX().contains(i));
				assertEquals(seqs.get(i).charAt(pos), cnode.getCH());
				pos++;
			}
		}
	}
	
}
