package jseqa.poa.cli;

import static org.junit.Assert.assertEquals;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.ParseException;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class POAlignAlignmentInputTest {

    @Rule
    public TemporaryFolder temporaryFolder = new TemporaryFolder();

    @Test(expected = ParseException.class)
    public void fastaAndAlignmentInputsAreMutuallyExclusive() throws Exception {
        new DefaultParser().parse(POAlign.createOptions(),
                new String[] {"--input", "reads.fasta", "--alignment", "reads.aln"});
    }

    @Test
    public void readsMultipleClustalBlocksInIdentifierOrder() throws Exception {
        Path input = write("blocks.aln",
                "CLUSTAL W multiple sequence alignment",
                "",
                "read1  AC-G  3",
                "read2  ACTG  4",
                "       ** *",
                "",
                "read1  TT",
                "read2  T-",
                "       * ");
        assertEquals(Arrays.asList("AC-GTT", "ACTGT-"),
                POAlign.readClustalAlignment(input));
    }

    @Test
    public void alignmentOptionImportsMsaWithoutRealigningIt() throws Exception {
        Path input = write("input.aln",
                "CLUSTAL W multiple sequence alignment",
                "",
                "reference  AC-GT",
                "variant    ACTGT");
        Path output = temporaryFolder.newFile("output.aln").toPath();
        POAlign.main(new String[] {
                "--alignment", input.toString(),
                "--type", "fitting",
                "--output", output.toString()
        });
        assertEquals(Arrays.asList("AC-GT", "ACTGT"),
                Files.readAllLines(output, StandardCharsets.UTF_8));
    }

    @Test(expected = IllegalArgumentException.class)
    public void rejectsFilesWithoutClustalHeader() throws Exception {
        Path input = write("not-clustal.aln", "read1 AC-GT", "read2 ACTGT");
        POAlign.readClustalAlignment(input);
    }

    private Path write(String filename, String... lines) throws Exception {
        Path path = temporaryFolder.newFile(filename).toPath();
        List<String> content = Arrays.asList(lines);
        Files.write(path, content, StandardCharsets.UTF_8);
        return path;
    }
}
