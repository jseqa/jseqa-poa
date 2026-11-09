package jseqa.poa.cli;

import com.sun.management.ThreadMXBean;
import java.lang.management.GarbageCollectorMXBean;
import java.lang.management.ManagementFactory;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.List;
import jseqa.fastq.Entry;
import jseqa.fastq.FastaReader;
import jseqa.poa.POGNode;
import jseqa.poa.POGraphLocal;
import jseqa.poa.PositionedPOGNodeFactory;

/** Repeatable single-thread POA benchmark excluding input and graph setup. */
public final class POAPerformanceBenchmark {
    private static final String DEFAULT_REFERENCE =
            "/Users/stef/git-self/igpoa/human-igh-anarci-refs/V.aln";
    private static final String DEFAULT_QUERIES =
            "/Users/stef/git-self/igpoa/atb_example_data/human_igh_50ents.fna";

    public static void main(String[] args) throws Exception {
        Path reference = Paths.get(args.length > 0 ? args[0] : DEFAULT_REFERENCE);
        Path queryFile = Paths.get(args.length > 1 ? args[1] : DEFAULT_QUERIES);
        int measuredPasses = args.length > 2 ? Integer.parseInt(args[2]) : 5;

        List<String> alignment = POAlign.readClustalAlignment(reference);
        List<String> queries = readFasta(queryFile);
        POGraphLocal graph = new POGraphLocal(alignment, new PositionedPOGNodeFactory());

        List<String> expectedPaths = alignAll(graph, queries);
        for (int warmup = 0; warmup < 3; warmup++) {
            assertSamePaths(expectedPaths, alignAll(graph, queries));
        }

        ThreadMXBean allocationBean = allocationBean();
        long threadId = Thread.currentThread().getId();
        long allocatedBefore = allocationBean == null ? -1
                : allocationBean.getThreadAllocatedBytes(threadId);
        long collectionsBefore = collectionCount();
        long started = System.nanoTime();

        for (int pass = 0; pass < measuredPasses; pass++) {
            assertSamePaths(expectedPaths, alignAll(graph, queries));
        }

        long elapsedNanos = System.nanoTime() - started;
        long collections = collectionCount() - collectionsBefore;
        long totalReads = (long) measuredPasses * queries.size();
        double elapsedSeconds = elapsedNanos / 1_000_000_000.0;

        System.out.printf("reference rows: %d, graph nodes: %d%n",
                alignment.size(), graph.numNodes());
        System.out.printf("queries: %d, measured passes: %d%n",
                queries.size(), measuredPasses);
        System.out.printf("output SHA-256: %s%n", digest(expectedPaths));
        System.out.printf("wall time: %.3f s%n", elapsedSeconds);
        System.out.printf("throughput: %.1f reads/s%n", totalReads / elapsedSeconds);
        if (allocationBean != null) {
            long allocated = allocationBean.getThreadAllocatedBytes(threadId) - allocatedBefore;
            System.out.printf("allocated: %.1f KiB/read%n", allocated / 1024.0 / totalReads);
        } else {
            System.out.println("allocated: unavailable on this JVM");
        }
        System.out.printf("garbage collections: %d%n", collections);
    }

    private static List<String> alignAll(POGraphLocal graph, List<String> queries) {
        List<String> paths = new ArrayList<String>(queries.size());
        for (String query : queries) {
            StringBuilder path = new StringBuilder();
            for (POGNode node : graph.alignToGraph(query)) {
                path.append(node == null ? '-' : node.getCH());
            }
            paths.add(path.toString() + "@" + graph.queryInterval());
        }
        return paths;
    }

    private static void assertSamePaths(List<String> expected, List<String> actual) {
        if (!expected.equals(actual)) {
            throw new AssertionError("alignment output changed during benchmark");
        }
    }

    private static String digest(List<String> paths) throws Exception {
        MessageDigest digest = MessageDigest.getInstance("SHA-256");
        for (String path : paths) {
            digest.update(path.getBytes(StandardCharsets.UTF_8));
            digest.update((byte) '\n');
        }
        StringBuilder hex = new StringBuilder();
        for (byte value : digest.digest()) {
            hex.append(String.format("%02x", value & 0xff));
        }
        return hex.toString();
    }

    private static List<String> readFasta(Path path) throws Exception {
        List<String> sequences = new ArrayList<String>();
        try (FastaReader reader = new FastaReader(path)) {
            for (Entry entry : reader) {
                sequences.add(entry.sequence());
            }
        }
        return sequences;
    }

    private static ThreadMXBean allocationBean() {
        java.lang.management.ThreadMXBean bean = ManagementFactory.getThreadMXBean();
        if (!(bean instanceof ThreadMXBean)) {
            return null;
        }
        ThreadMXBean allocationBean = (ThreadMXBean) bean;
        if (!allocationBean.isThreadAllocatedMemorySupported()) {
            return null;
        }
        if (!allocationBean.isThreadAllocatedMemoryEnabled()) {
            allocationBean.setThreadAllocatedMemoryEnabled(true);
        }
        return allocationBean;
    }

    private static long collectionCount() {
        long count = 0;
        for (GarbageCollectorMXBean bean : ManagementFactory.getGarbageCollectorMXBeans()) {
            if (bean.getCollectionCount() >= 0) {
                count += bean.getCollectionCount();
            }
        }
        return count;
    }

    private POAPerformanceBenchmark() {}
}
