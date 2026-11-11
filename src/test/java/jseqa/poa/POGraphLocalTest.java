package jseqa.poa;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.junit.Test;

import jseqa.alignment.substitution.CustomDnaSubstitutionMatrix;
import jseqa.util.Interval;


public class POGraphLocalTest {

	@Test
	public void testSimple() {
		POGraphLocal pog = new POGraphLocal("ACGCTT");

		pog.add("ACCCTT");
		pog.add("ACTCTT");
		pog.add("ACCTT");
		pog.add("AACCCTTTT");
		pog.computeInternalDS();
		
		String cons = pog.consensus();
//		System.out.println(cons);
//		pog.toDOT("/tmp/simple_graph.dot");

		String exp = "ACCCTT";
		assertEquals(exp, cons);
		
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
		//  AAGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGGTTTGTGGAGGCAA
		String exp = "AAGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGGTTTGTGGAGGCAA";
		POGraphLocal pog = new POGraphLocal();

		for(int i = 0; i < seqs.size(); i++) {
//			System.out.println(i + "\t" + seqs.get(i));
			pog.add(seqs.get(i));
		}
		pog.computeInternalDS();
//		pog.toDOT("test_graph.dot");
		String cons = pog.consensus();
		List<String> algn = pog.getAlignment();
//		for(String s : algn) {
//			System.out.println(s);
//		}
//		assertEquals(exp, cons);
		
		int num = 0;
		for(POGNode currNode : pog) {
			num++;
		}
		// The final graph/query boundary is a valid alignment column.
		assertEquals(49, num);
		
		for(int i = 0; i < seqs.size(); i++) {
			List<POGNode> path = pog.alignToGraph(seqs.get(i));
			String substr = path.stream().map(n -> ""+n.getCH()).collect(Collectors.joining(""));			
			int substrIdx = seqs.get(i).indexOf(substr);
			assertTrue(substrIdx >= 0);
		}
	}
	
	@Test
	public void testAlignWithoutAdding() {
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

//		POGraph pog = new POGraph(); //seqs.get(0));
		POGraphGlobal pog = new POGraphGlobal(); //seqs.get(0));
//		for(int i = 0; i < seqs.size(); i++) {
		for(int i = 0; i < 3; i++) {
			pog.add(seqs.get(i));
		}
//		pog.toDOT("/tmp/small_aln_pre.dot");
		pog.computeInternalDS();
//		pog.toDOT("/tmp/small_aln.dot");
		String cons = pog.consensus();
		List<String> algn = pog.getAlignment();
		
		String s1 = seqs.get(0);
		pog.alignSequence(s1);
		List<POGNode> p1 = pog.backtrackPath(s1);
		for(int i = 0; i < p1.size(); i++) {
//			System.out.println(p1.get(i));
			assertEquals(s1.charAt(i), p1.get(i).getCH());
		}

		// make changes and check for null entries in path
		StringBuilder s1Mod = new StringBuilder(s1);
		s1Mod.setCharAt(1, 'A');
		s1Mod.setCharAt(45, 'T');
		pog.alignSequence(s1Mod.toString());
		List<POGNode> p2 = pog.backtrackPath(s1Mod.toString());
		for(int i = 0; i < p2.size(); i++) {
			if(i != 1 && i != 45) {
				assertEquals(s1.charAt(i), p2.get(i).getCH());
			} else {
				assertEquals(null, p2.get(i));
			}
		}
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
			add("IGHJ5*02");
			add("IGHJ6*01");
			add("IGHJ6*02");
		}};
		
		CustomDnaSubstitutionMatrix scoremat = new CustomDnaSubstitutionMatrix(5, -4, -4, -1);
//		CustomDNAScoreMatrix scoremat = new CustomDNAScoreMatrix(5, -20, -20, -1);
//		CustomDNAScoreMatrix scoremat = new CustomDNAScoreMatrix(5, -4, -10, -1);
//		POGraphLocal poa = new POGraphLocal(scoremat);
//		POGraph poaG = new POGraph(scoremat);
		POGraphFitting poaG = new POGraphFitting(scoremat);
//		POGraph poa = new POGraph(scoremat);
		
//		poa.addAll(refSeqs);
		for(String s : refSeqs) {
			poaG.add(s);
//			System.out.println(poaG);
		}
		
		POGraphLocal poa = new POGraphLocal(poaG);
//		POGraph poa = poaG;

		int i = 0;
		for(String jseq : refSeqs) {
			List<POGNode> jpath = poa.alignToGraph(jseq);
			for(POGNode node : jpath) {
				assertTrue(node.getIDX().contains(i));
			}
			i++;
		}
		
		
		String seq = "AATAGCTACGATGACTATGGTGATTACATGACTGGTTGGATCTCTGGGGCCAGGGCACCCTGGTCACCGTCTCTTCAN";
		List<POGNode> path = poa.alignToGraph(seq);		
//		poa.toDOT("/tmp/tmpgraph.dot");
		Interval qint = poa.queryInterval();
		
		String subSeq = seq.substring(qint.start(), qint.end());
//		System.out.println(seq);
		String buf = StringUtils.repeat(' ', qint.start());
//		System.out.println(buf + subSeq);
//		System.out.println(path);
//		System.out.println(path.size());
		Interval jint = new Interval(seq.indexOf(refSeqs.get(7)), seq.indexOf(refSeqs.get(7))+refSeqs.get(7).length());
//		System.out.println(qint + "\n" + jint);
	}	
	
	@Test
	public void testVExample() {
		List<String> refSeqs = new ArrayList<String>() {{
			add("CAGGAGCAGCTGAAGGAGTCCGGGGGTCGTCTGGTCATGCCTGGAGGAATCCTGACACTCACCTGCACAGCCTCTGGATCCAACATCAGTAGCTATGGAGTGAGCTGGTTCCGCCAGGCTCCAGGGAAGGGGCTCGAATGGATAAGATACATTAGTTATGGTGGTAGCGCATACTACAAGAGCTGGGTGAAAGGCCGATTCACCATCTCCAAAACCTCGAGCACGGTGGATCTGAAAATGACCAGTCTGACAGCCTCGGACACGGCCACCTATTTCTGTGCGAGAGA");
			add("CAGGAGCAGCTGGAGGAATCCGGGGGAGGCCTGGTCACGCCTGGAGGAACCCTGACACTCACCTGCACAGTCTCTGGATTCTCCCTCAGTAGCTATGGAGTAAGCTGGGTCCGCCAGGCTGCAGGGAAGGGGCTGGAATGGATCGGATACATTAGTAGTAGTGGTAGCGCATACTACGCGAGTTGGGTAAATGGTCGATTCACCATCTCCAAAACCTCGAGCACGGTGGATCTGAAAATGACCAGTCTGAGAGCCGCGGACACGGCCACCTATTTCTGTGCCAGAAA");
			add("CAGGAGCAGCAGAAGGAGTCCGGGGGTCGACTGGTCATGCCTGGAGGATCCCTGACACTCACCTGCACAGTCTCTGGATTCTCCCTCAGTAGCTACAACATGGGCTGGGTCCGCCAGGCTCCAGGGGAGGGGCTGGAATACATCGGATGGATTAGTACTGGTGGTAGCGCATACTACGCGAGCTGGGTGAATGGTCGATTCACCATCTCCAAAACCTCGACCACGATGGATCTGAAAATGACCAGTCTGACAGCCGCGGACACGGCCACCTATTTCTGTGCCAGAGA");
			add("GAGGATCAGCTAGTGGAGTCCAGGGGAGGCCTGCTCCAGCCTGGGGGATCCCTGACACTCTCCTGCAAAGCCTCTGGATTTACCATCAGTAGCTATGGAGTGAGCTGGGTCTGCCAGGCTCCAGGGAAGGGGCTGGAGTGGATCGGATACATTTATACTGATGGTAGCACATACTACGCGAGCTGGGTGAATGGCCGATTCACCATCTCCAAAACCTCGACCACGGTGGATCTGAAAATGACCAGTCTGACAGCCGCGGACACGGCCACCTATTTCTGTGCAAGAAA");
			add("CAGGAGCAGCTGAAGGAGTCCGGGGGAGGCCTGGTCACGCCTGGAGGAATCCTGTCACTCACCTGCACAGCCTCTGGATTCTCCATCAGTAGCTACAGAATGGGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAATACATCGGATACATTAGTTACGGTGGTAGCGCATACTACAAGAGCTGGGTGAAAGGCCGATTCACCATCTCCAAAACCTCAAGCACGGTGGATCTGAAAATGACCAGTCTGACAGCCTCGGACAAGGCCACCTATTTCTGTGCCAGAGA");
			add("CAGGAGCAGCTGAAGGAGTCCGGGGGAGTCCTGGTCACGCCTGGAGGAATCCTGTCACTCACCTGCACAGCCTCTGGATTCTCCATCAGTAGCTACAGAATGGGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAATACATCGGAATCATTTATACTGGTGGTAGCGCATACTACGCGAGCTGGGTGAATGGTCGATTCACCATCTCCAAAACCTCGAGCACGGTGGATCTGAAAATGACCAGTCTGACAGCCGCGGACATGGCCACCTATTTCTGTGCCAGAGA");
			add("CAGTCAGTGAAGGAGTCCGAGGGAGGTCTCTTCAAGCCAACGGATACCCTGACACTCACCTGCACAGTCTCTGGATTCTCCCTCAGTAGCTATGCAATAAGCTGGGTCCGCCAGGCTCCAGGGAACGGGCTGGAATGGATCGGAATCATTAATAGTTATGGTAGCACATACTACGCGAGCTGGGCGAAAAGCCGATCCACCATCACCAGAAACACCAACGAGAACACGGTGACTCTGAAAATGACCAGTCTGACAGCCGCGGACACGGCCACCTATTTCTGTGCGAGAGA");
			add("CAGTCGCTGGAGGAGTCCAGGGGAGGCCTGATCAAGCCTGGAGGAACCCTGACACTCACCTGCACAGCCTCTGGATTCACCATCAGTAGCTACGACATGAGCTGGGTCCGACAGGCTCCAGGGAAGGAGCTGGAGTGGATCGGATACATTAGTTATGGTGGTAGCGCATACTACGCGAGCTGGGCGAAAAGCCGATCCACCATCACCAGAAACACCAACGAGAACACGGTGACTCTGAAAATGACCAGTCTGACAGCCGCGGACACGGCCACCTATTTCTGTGCGAGAGA");
			add("CAGGAGCAGCTGAAGGAGTCCGGAAGAGGCCTGGTCCAGCCTGGGGGATCCCTGACACTCACCTGCACAGCCTCCGGATTCACCATCAGTAGCAATGCAATGATCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAATACATCGGATTCATTAATACTGGTGGTAGCACTAACTACGCGAGCTGGGCGAATGGCCGATTCACCATCTCCAGCGACAACGCCCAGAACACGGTGTCTCTGCAAATGAACAGTCTGACAGCCGCGGACACCTATTTCTGTGCCAGAGA");
		}};
		List<String> ids = new ArrayList<String>() {{
			add("IGHV1S1*01");
			add("IGHV1S13*01");
			add("IGHV1S17*01");
			add("IGHV1S21*01");
			add("IGHV1S24*01");
			add("IGHV1S25*01");
			add("IGHV1S26*01");
			add("IGHV1S28*01");
			add("IGHV1S29*01");
		}};
		CustomDnaSubstitutionMatrix scoremat = new CustomDnaSubstitutionMatrix(5, -4, -4, -1);
//		POGraph poaG = new POGraph(scoremat);
		POGraphLocal poaG = new POGraphLocal(scoremat);
		for(int i = 0; i < ids.size(); i++) {
//		for(int i = 0; i < 2; i++) {
//			poaG.add(refSeqs.get(i).substring(100));
			poaG.add(refSeqs.get(i));
		}
//		poaG.toDOT("/tmp/vsmall_pre.dot");
		poaG.computeInternalDS();
//		poaG.toDOT("/tmp/vsmall_post.dot");
		
		
		List<POGNode> path = poaG.alignToGraph(refSeqs.get(8).substring(100));
		for(POGNode node : path) {
			boolean hasRead9 = (node != null ? node.getIDX().contains(8) : false);
//			System.out.println(hasRead9 + "\t" + node);
			assertTrue(hasRead9);
		}
//		System.out.println(refSeqs.get(8).substring(100));
	}
	
}
