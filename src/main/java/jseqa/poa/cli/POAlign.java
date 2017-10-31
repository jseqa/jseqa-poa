package jseqa.poa.cli;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import jseqa.alignment.substitution.CustomDnaSubstitutionMatrix;
import jseqa.fastq.Entry;
import jseqa.fastq.FastaReader;
import jseqa.poa.POGraphFitting;
import jseqa.poa.POGraphGlobal;
import jseqa.poa.POGraphLocal;
import jseqa.poa.POGraphOverlap;
import jseqa.poa.PartialOrderAlignment;


public class POAlign {
	
	private final static Logger log = LogManager.getLogger(POAlign.class);
	
	public static void main(String[] args) {		
		CommandLineParser clParser = new DefaultParser();		
		Options opts = new Options();
		opts.addOption("h", "help", false, "display help");
		opts.addOption("v", "verbose", false, "verbose output");
		opts.addOption("i", "input", true, "input FASTA file");
		opts.addOption("o", "output", true, "output MSA file");
		opts.addOption("d", "dot", true, "output to DOT file");
		opts.addOption("t", "type", true, "type of alignment [global/local/fitting/overlap]");
		HelpFormatter formatter = new HelpFormatter();
		Path inputFasta = null;
		Path outFile = null;
		Path outDotFile = null;
		String alnType = "global"; 
		try {
			CommandLine cmd = clParser.parse(opts, args);
			
			if( cmd.hasOption("help") || (cmd.getOptions().length == 0)) {
				formatter.printHelp("help", opts);
				System.exit(-1);
			}
			if( cmd.hasOption("1") && !cmd.hasOption("2") ) {
				System.err.println("Missing second read file, use -2");
				System.exit(-1);
			}
			if(cmd.hasOption("i")) {
				inputFasta = Paths.get(cmd.getOptionValue("i"));
			} else {
				System.err.println("input file required (-i)");
				System.exit(-1);
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
					alnType = "global";
				}
			}
			if(cmd.hasOption("o")) {
				outFile = Paths.get(cmd.getOptionValue("o"));
			}
			if(cmd.hasOption("d")) {
				outDotFile = Paths.get(cmd.getOptionValue("d"));
			}
		} catch (ParseException e) {
			e.printStackTrace();
		}
		CustomDnaSubstitutionMatrix scoremat = new CustomDnaSubstitutionMatrix(5, -4, -4, -1);
		PartialOrderAlignment poa = null;
		
		if(alnType.equalsIgnoreCase("global")) {
			poa = new POGraphGlobal(scoremat);
		} else if(alnType.equalsIgnoreCase("local")) {
			poa = new POGraphLocal(scoremat);
		} else if(alnType.equalsIgnoreCase("fitting")) {
			poa = new POGraphFitting(scoremat);
		} else if(alnType.equalsIgnoreCase("overlap")) {
			poa = new POGraphOverlap(scoremat);
		} else {
			poa = new POGraphGlobal(scoremat);
		}
		
		
		try(FastaReader faread= new FastaReader(inputFasta)) {
			for(Entry ent : faread) {
				poa.add(ent.sequence());
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		List<String> alnRows = poa.getAlignment();
		if(outFile == null) {
			System.out.println(alnRows);
		} else {
			System.err.println("WRITING TO FILE NOT YET SUPPORTED");
		}
		
		if(outDotFile != null) {
			poa.toDOT(outDotFile.toString());
		}
	}
	
	private POAlign() {}
}

