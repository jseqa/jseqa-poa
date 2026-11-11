package jseqa.poa.cli;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.net.URL;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import jseqa.poa.POGraphFitting;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class POAlignOutputTest {
	private static final Pattern DOT_NODE = Pattern.compile(
			"^(\\d+) \\[shape=circle,.*label=\\\"(.)\\\".*\\];$");
	private static final Pattern DOT_EDGE = Pattern.compile("^(\\d+) -> (\\d+);$");

    @Rule
    public TemporaryFolder temporaryFolder = new TemporaryFolder();

    @Test
    public void writesOneAlignmentRow() throws Exception {
        Path output = temporaryFolder.newFile("one.aln").toPath();
        POAlign.writeAlignment(output, Collections.singletonList("ACGT"));
        assertEquals(Collections.singletonList("ACGT"), Files.readAllLines(output, StandardCharsets.UTF_8));
    }

    @Test
    public void writesEachAlignmentRowOnSeparateLine() throws Exception {
        Path output = temporaryFolder.newFile("many.aln").toPath();
        POAlign.writeAlignment(output, Arrays.asList("AC-G", "A-TG", "ACTG"));
        assertEquals(Arrays.asList("AC-G", "A-TG", "ACTG"),
                Files.readAllLines(output, StandardCharsets.UTF_8));
    }

    @Test
    public void replacesExistingOutputInsteadOfAppending() throws Exception {
        Path output = temporaryFolder.newFile("replace.aln").toPath();
        Files.write(output, Collections.singletonList("old"), StandardCharsets.UTF_8);
        POAlign.writeAlignment(output, Collections.singletonList("new"));
        assertEquals(Collections.singletonList("new"), Files.readAllLines(output, StandardCharsets.UTF_8));
    }

	@Test
	public void buildsPoaFromExpectedAlignmentAndWritesCorrectDotGraph() throws Exception {
		URL alignmentResource = POAlignOutputTest.class.getResource("cli-test.expected.aln");
		assertNotNull("missing cli-test.expected.aln test resource", alignmentResource);
		Path alignment = Paths.get(alignmentResource.toURI());
		List<String> rows = Files.readAllLines(alignment, StandardCharsets.UTF_8);
		POGraphFitting graph = new POGraphFitting(rows);
		Path dotOutput = temporaryFolder.newFile("cli-test.dot").toPath();

		graph.toDOT(dotOutput.toString(), true);

		List<String> dot = Files.readAllLines(dotOutput, StandardCharsets.UTF_8);
		assertEquals("digraph g {", dot.get(0));
		assertEquals("}", dot.get(dot.size() - 1));

		Set<Integer> nodeIds = new HashSet<Integer>();
		Map<Character,Integer> residueCounts = new HashMap<Character,Integer>();
		int edgeCount = 0;
		for (String line : dot.subList(1, dot.size() - 1)) {
			Matcher node = DOT_NODE.matcher(line);
			Matcher edge = DOT_EDGE.matcher(line);
			if (node.matches()) {
				nodeIds.add(Integer.valueOf(node.group(1)));
				char residue = node.group(2).charAt(0);
				residueCounts.put(residue, residueCounts.getOrDefault(residue, 0) + 1);
			} else if (edge.matches()) {
				edgeCount++;
				assertTrue("edge source must be declared", nodeIds.contains(Integer.valueOf(edge.group(1))));
				assertTrue("edge target must be declared", nodeIds.contains(Integer.valueOf(edge.group(2))));
			} else {
				throw new AssertionError("Unexpected DOT line: " + line);
			}
		}
		assertEquals(10, graph.numNodes());
		assertEquals(graph.numNodes(), nodeIds.size());
		assertEquals(graph.getGraph().edgeSet().size(), edgeCount);
		assertEquals(Integer.valueOf(2), residueCounts.get('A'));
		assertEquals(Integer.valueOf(2), residueCounts.get('C'));
		assertEquals(Integer.valueOf(3), residueCounts.get('G'));
		assertEquals(Integer.valueOf(3), residueCounts.get('T'));
	}
}
