package jseqa.poa;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Optional;

import org.junit.Test;

import jseqa.util.Pair;


public class POGPtrTest {

	@Test
	public void testPOGPtr() {
		POGPtr p1 = new POGPtr(new Pair<Integer,Integer>(1,1), 5);
		POGPtr p2 = new POGPtr(new Pair<Integer,Integer>(1,2), 10);
		POGPtr p3 = new POGPtr(new Pair<Integer,Integer>(1,3), -5);
		
		ArrayList<POGPtr> lst = new ArrayList<POGPtr>();
		lst.add(p1);
		lst.add(p2);
		lst.add(p3);
		
		Optional<POGPtr> curr = lst.stream().max(Comparator.comparing(i -> i.score()));
		assertEquals(10, curr.get().score());
	}

}
