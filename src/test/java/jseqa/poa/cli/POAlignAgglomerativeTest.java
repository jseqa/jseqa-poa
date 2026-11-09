package jseqa.poa.cli;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.cli.Options;
import org.junit.Test;

public class POAlignAgglomerativeTest {

    @Test
    public void commandLineExposesAgglomerativeOptions() {
        Options options = POAlign.createOptions();
        assertTrue(options.hasOption("progressive"));
        assertTrue(options.hasOption("linkage"));
        assertTrue(options.hasOption("pairwise-mode"));
        assertTrue(options.hasOption("score-normalization"));
        assertTrue(options.hasOption("guide-tree-output"));
        assertTrue(options.hasOption("threads"));
    }

    @Test
    public void agglomerativeCliWritesPreservedMsaAndGuideTree() throws Exception {
        Path fasta = Files.createTempFile("poa-agglomerative", ".fasta");
        Path alignment = Files.createTempFile("poa-agglomerative", ".aln");
        Path tree = Files.createTempFile("poa-agglomerative", ".nwk");
        List<String> sequences = Arrays.asList("TGCGAGA", "CAACAGA", "TGCGAGA");
        Files.write(fasta, Arrays.asList(
                ">first", sequences.get(0),
                ">second", sequences.get(1),
                ">duplicate", sequences.get(2)), StandardCharsets.UTF_8);

        POAlign.main(new String[] {
                "--input", fasta.toString(),
                "--output", alignment.toString(),
                "--progressive", "agglomerative",
                "--linkage", "upgma",
                "--pairwise-mode", "global",
                "--score-normalization", "self-score",
                "--threads", "2",
                "--guide-tree-output", tree.toString()
        });

        List<String> rows = Files.readAllLines(alignment, StandardCharsets.UTF_8);
        assertEquals(3, rows.size());
        List<String> reconstructed = ungapped(rows);
        assertEquals(2, java.util.Collections.frequency(reconstructed, "TGCGAGA"));
        assertEquals(1, java.util.Collections.frequency(reconstructed, "CAACAGA"));
        String newick = new String(Files.readAllBytes(tree), StandardCharsets.UTF_8).trim();
        assertTrue(newick.endsWith(";"));
        assertTrue(newick.contains("0"));
        assertTrue(newick.contains("1"));
        assertTrue(newick.contains("2"));
    }

    @Test
    public void inputOrderRemainsTheDefaultProgressiveMode() throws Exception {
        Path fasta = Files.createTempFile("poa-input-order", ".fasta");
        Path alignment = Files.createTempFile("poa-input-order", ".aln");
        Files.write(fasta, Arrays.asList(">one", "ACGT", ">two", "AGGT"),
                StandardCharsets.UTF_8);
        POAlign.main(new String[] {
                "--input", fasta.toString(), "--output", alignment.toString()
        });
        assertEquals(Arrays.asList("ACGT", "AGGT"),
                ungapped(Files.readAllLines(alignment, StandardCharsets.UTF_8)));
    }

    private static List<String> ungapped(List<String> rows) {
        java.util.ArrayList<String> sequences = new java.util.ArrayList<String>();
        for (String row : rows) sequences.add(row.replace("-", ""));
        return sequences;
    }
}
