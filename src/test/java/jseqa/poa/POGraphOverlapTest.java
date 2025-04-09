package jseqa.poa;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.junit.Test;

import jseqa.alignment.substitution.CustomDnaSubstitutionMatrix;
import jseqa.util.Interval;


public class POGraphOverlapTest {
	
	@Test
	public void testSimple() {
		POGraphOverlap pog = new POGraphOverlap();

		List<String> seqs = Arrays.asList(
		"AACCCTTTT",
		"ACGCTT",
		"ACCCTT",
		"ACTCTT",
		"ACCTT");
		for(String seq : seqs) {
			pog.add(seq);
		}
		
		String cons = pog.consensus();
		String exp = "AACCCTTTT";
		assertEquals(exp, cons);
		
		for(String inseq : seqs) {
			List<POGNode> path = pog.alignToGraph(inseq);
			for(int i = 0; i < inseq.length(); i++) {
				assertEquals(inseq.charAt(i), path.get(i).getCH());
			}
		}
				
		int num = 0;
		for(POGNode currNode : pog) {
			num++;
		}
		assertEquals(exp.length(), num);
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
		
		
		POGraphOverlap pog = new POGraphOverlap(new CustomDnaSubstitutionMatrix(5, -4, -6, -4));
		
		for(int i = 0; i < seqs.size(); i++) {
			pog.add(seqs.get(i));
		}
		
		// should retrieve full 4th seq
		List<POGNode> path4 = pog.alignToGraph(seqs.get(4));
		for(int i = 0; i < seqs.get(4).length(); i++) {
			assertEquals(seqs.get(4).charAt(i), path4.get(i).getCH());
		}
		
		String query = "CTGCGTTTGCTGGACTGATGTACTTGT";		
		List<POGNode> qPath = pog.alignToGraph(query);
		for(int i = 0; i < query.length(); i++) {
			assertEquals(query.charAt(i), qPath.get(i).getCH());
		}
		
		String cons = pog.consensus();
		List<String> algn = pog.getAlignment();		
		String exp = "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGGTTTGTGGAGGCAA";
		assertEquals(exp, cons);
		
		int num = 0;
		for(POGNode currNode : pog) {
			num++;
		}
		assertEquals(48, num);
		//pog.toDOT("test-align-pog.dot", true);
	}

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
			add("ATGACTGGTTGGATCTCTGGGGCCAGGGCACCCTGGTCACCGTCTCTTCAN");		//IGHJ5*01
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
		
		CustomDnaSubstitutionMatrix scoremat = new CustomDnaSubstitutionMatrix(5, -8, -4, -1);
//		POGraphOverlap poa = new POGraphOverlap(scoremat);
//		POGraph poa = new POGraph(scoremat);
		POGraphGlobal poa = new POGraphGlobal(scoremat);

//		poa.addAll(refSeqs);
		for(String s : refSeqs) {
			poa.add(s);
		}

//		System.out.println(poa);
		
		POGraphOverlap poaO = new POGraphOverlap(poa);
		PartialOrderAlignment poaI = poaO;
		
		String seq = "GAGCATATACTAGTAGTAGTGGTTATTATATACTGATGCTTTTGATCCCTGGGGCCCAGGCACCCTGGTCACCGTCTCCTCAG";
//		String seq = "AGTTACTATAGTAGTGGCTGGGGTGATGACTGGTTGGATCTCTGGGGCCAGGGCACCCTGGTCACCGTCTCTTCAN";
		
		List<POGNode> path = poaI.alignToGraph(seq);
//		poaI.toDOT("/tmp/tmpgraph.dot");
		Interval qint = poaI.queryInterval();
		String subSeq = seq.substring(qint.start(), qint.end());
//		System.out.println(seq);
		String buf = StringUtils.repeat(' ', qint.start());
//		System.out.println(buf + subSeq);
//		System.out.println(path);
//		System.out.println(path.size());
		Interval jint = new Interval(seq.indexOf(refSeqs.get(7)), seq.indexOf(refSeqs.get(7))+refSeqs.get(7).length());
//		System.out.println(qint + "\n" + jint);
	}
}
