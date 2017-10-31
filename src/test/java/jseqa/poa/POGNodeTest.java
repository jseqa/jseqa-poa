package jseqa.poa;

import static org.junit.Assert.*;

import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

public class POGNodeTest {

	@Test
	public void testPOGNodeCharInt() {
		POGNode p1 = new POGNode('A', 0);
		POGNode p2 = new POGNode('A', 0);
		POGNode p3 = new POGNode('A', 0);
		
		
		assertEquals(1, p2.getCounter()-p1.getCounter());
		assertEquals(2, p3.getCounter()-p1.getCounter());
	}

	@Test
	public void testHash() {
		POGNode p1 = new POGNode('A', 0);
		POGNode p2 = new POGNode('A', 0);
		POGNode p3 = new POGNode('A', 0);
		
		Set<POGNode> h = new HashSet<POGNode>();
		h.add(p1);
		h.add(p2);
		h.add(p3);
		
		assertEquals(3, h.size());
	}
}
