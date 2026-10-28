package jseqa.poa.cli;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import jseqa.alignment.substitution.CustomDnaSubstitutionMatrix;
import jseqa.fastq.Entry;
import jseqa.fastq.FastaReader;
import jseqa.poa.POGraphFitting;
import jseqa.poa.POGraphGlobal;
import jseqa.poa.POGraphLocal;
import jseqa.poa.POGraphOverlap;
import jseqa.poa.PartialOrderAlignment;


public class POAlign {
	
	public static void main(String[] args) {		
		CommandLineParser clParser = new DefaultParser();		
		Options opts = createOptions();
		HelpFormatter formatter = new HelpFormatter();
		Path inputFasta = null;
		Path inputAlignment = null;
		Path outFile = null;
		Path outDotFile = null;
		String alnType = "global"; 
		try {
			CommandLine cmd = clParser.parse(opts, args);
			
			if( cmd.hasOption("help") || (cmd.getOptions().length == 0)) {
				formatter.printHelp("poalign", opts);
				return;
			}
			if(cmd.hasOption("i")) {
				inputFasta = Paths.get(cmd.getOptionValue("i"));
			} else if (cmd.hasOption("a")) {
				inputAlignment = Paths.get(cmd.getOptionValue("a"));
			} else {
				throw new ParseException("one input is required: --input or --alignment");
			}
			//
			if(cmd.hasOption("t")) {
				String val = cmd.getOptionValue("t");
				if(val.equalsIgnoreCase("global")) {
					alnType = "global";
				} else if(val.equalsIgnoreCase("local")) {
					alnType = "local";
				} else if(val.equalsIgnoreCase("fitting")) {
					alnType = "fitting";
				} else if(val.equalsIgnoreCase("overlap")) {
					alnType = "overlap";
				} else {
					throw new ParseException("Unknown alignment type: " + val);
				}
			}
			if(cmd.hasOption("o")) {
				outFile = Paths.get(cmd.getOptionValue("o"));
			}
			if(cmd.hasOption("d")) {
				outDotFile = Paths.get(cmd.getOptionValue("d"));
			}
		} catch (ParseException e) {
			System.err.println(e.getMessage());
			formatter.printHelp("poalign", opts);
			return;
		}
		CustomDnaSubstitutionMatrix scoremat = new CustomDnaSubstitutionMatrix(5, -4, -4, -1);
		PartialOrderAlignment poa;

		try {
			if (inputAlignment != null) {
				poa = createGraphFromAlignment(alnType, scoremat,
						readClustalAlignment(inputAlignment));
			} else {
				poa = createEmptyGraph(alnType, scoremat);
				try(FastaReader faread = new FastaReader(inputFasta)) {
					for(Entry ent : faread) {
						poa.add(ent.sequence());
					}
				}
			}
		} catch (IOException | IllegalArgumentException e) {
			System.err.println("Unable to read input: " + e.getMessage());
			return;
		}
		
		List<String> alnRows = poa.getAlignment();
		if(outFile == null) {
			for (String row : alnRows) {
				System.out.println(row);
			}
		} else {
			try {
				writeAlignment(outFile, alnRows);
			} catch (IOException e) {
				System.err.println("Unable to write alignment: " + e.getMessage());
				return;
			}
		}
		
		if(outDotFile != null) {
			poa.toDOT(outDotFile.toString());
		}
	}

	static Options createOptions() {
		Options options = new Options();
		options.addOption("h", "help", false, "display help");
		options.addOption("v", "verbose", false, "verbose output");
		OptionGroup inputs = new OptionGroup();
		inputs.addOption(new Option("i", "input", true, "input unaligned FASTA file"));
		inputs.addOption(new Option("a", "alignment", true, "input CLUSTAL .aln file"));
		options.addOptionGroup(inputs);
		options.addOption("o", "output", true, "output MSA file");
		options.addOption("d", "dot", true, "output to DOT file");
		options.addOption("t", "type", true,
				"type of alignment [global/local/fitting/overlap]");
		return options;
	}

	static List<String> readClustalAlignment(Path input) throws IOException {
		List<String> lines = Files.readAllLines(input, StandardCharsets.UTF_8);
		Map<String,StringBuilder> rowsById = new LinkedHashMap<String,StringBuilder>();
		boolean foundHeader = false;
		for (String line : lines) {
			String trimmed = line.trim();
			if (trimmed.isEmpty()) {
				continue;
			}
			if (!foundHeader) {
				String upper = trimmed.toUpperCase();
				if (!upper.startsWith("CLUSTAL") && !upper.startsWith("MUSCLE")) {
					throw new IllegalArgumentException("missing CLUSTAL header");
				}
				foundHeader = true;
				continue;
			}
			if (Character.isWhitespace(line.charAt(0)) || trimmed.startsWith("#")) {
				continue;
			}
			String[] fields = trimmed.split("\\s+");
			if (fields.length < 2 || !fields[1].matches("[A-Za-z*?.-]+")) {
				throw new IllegalArgumentException("invalid alignment line: " + line);
			}
			StringBuilder row = rowsById.get(fields[0]);
			if (row == null) {
				row = new StringBuilder();
				rowsById.put(fields[0], row);
			}
			row.append(fields[1].replace('.', '-'));
		}
		if (!foundHeader || rowsById.isEmpty()) {
			throw new IllegalArgumentException("alignment contains no sequence rows");
		}
		List<String> rows = new ArrayList<String>(rowsById.size());
		for (StringBuilder row : rowsById.values()) {
			rows.add(row.toString());
		}
		return rows;
	}

	private static PartialOrderAlignment createEmptyGraph(String alignmentType,
			CustomDnaSubstitutionMatrix scoremat) {
		if(alignmentType.equalsIgnoreCase("local")) {
			return new POGraphLocal(scoremat);
		} else if(alignmentType.equalsIgnoreCase("fitting")) {
			return new POGraphFitting(scoremat);
		} else if(alignmentType.equalsIgnoreCase("overlap")) {
			return new POGraphOverlap(scoremat);
		}
		return new POGraphGlobal(scoremat);
	}

	private static PartialOrderAlignment createGraphFromAlignment(String alignmentType,
			CustomDnaSubstitutionMatrix scoremat, List<String> rows) {
		if(alignmentType.equalsIgnoreCase("local")) {
			return new POGraphLocal(rows, scoremat);
		} else if(alignmentType.equalsIgnoreCase("fitting")) {
			return new POGraphFitting(rows, scoremat);
		} else if(alignmentType.equalsIgnoreCase("overlap")) {
			return new POGraphOverlap(rows, scoremat);
		}
		return new POGraphGlobal(rows, scoremat);
	}

	static void writeAlignment(Path output, List<String> rows) throws IOException {
		Files.write(output, rows, StandardCharsets.UTF_8);
	}
	
	private POAlign() {}
}

