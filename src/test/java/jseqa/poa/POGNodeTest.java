package jseqa.poa;

import static org.junit.Assert.*;

import java.util.HashSet;
import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

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

	@Test
	public void concurrentNodesHaveUniqueIdentifiers() throws Exception {
		final int threads = 8;
		final int nodesPerThread = 2000;
		final Set<Integer> identifiers = ConcurrentHashMap.newKeySet();
		ExecutorService executor = Executors.newFixedThreadPool(threads);
		List<Future<Boolean>> results = new ArrayList<Future<Boolean>>(threads);
		try {
			for (int thread = 0; thread < threads; thread++) {
				results.add(executor.submit(new Callable<Boolean>() {
					@Override
					public Boolean call() {
						for (int index = 0; index < nodesPerThread; index++) {
							if (!identifiers.add(new POGNode('A', 0).getCounter())) {
								return false;
							}
						}
						return true;
					}
				}));
			}
			for (Future<Boolean> result : results) {
				assertTrue(result.get());
			}
		} finally {
			executor.shutdownNow();
		}
		assertEquals(threads * nodesPerThread, identifiers.size());
	}
}
