package de.wuttke.vcf;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;

public class AddVcfDosage {    
	
	public static void main(String[] args) 
	throws Exception {
		AddVcfDosage tool = new AddVcfDosage();
		tool.execute(args);
	}

	private void execute(String[] args) throws IOException {
		parseCommandLine(args);
		prepareReadVcf();
		prepareWriteVcf();
        processVariants();
		cleanUp();
	}

	private void cleanUp() throws IOException {
		reader.close();
		writer.close();
	}

	private void processVariants() throws IOException {
		for (VariantContext variant : reader.iterator()) {
            progressLogger.record(variant.getContig(), variant.getStart());
            List<Genotype> gts = processVariantSamples(variant);
            writeVariant(variant, gts);
		}
	}

	private void writeVariant(VariantContext variant, List<Genotype> gts) {
		VariantContextBuilder variantBuilder = new VariantContextBuilder(variant);
		variantBuilder.genotypes(gts);
		VariantContext variant2 = variantBuilder.make();
		writer.add(variant2);
	}

	private List<Genotype> processVariantSamples(VariantContext variant) {
        List<Genotype> gts = new LinkedList<Genotype>();
		for (int sampleIdx = 0; sampleIdx < variant.getNSamples(); sampleIdx++) {
			Genotype g = variant.getGenotype(sampleIdx);
			processGenotype(variant, sampleIdx, g, gts);
		}
		
		return gts;
	}

	private void processGenotype(VariantContext variant, int sampleIdx, Genotype genotype, List<Genotype> gts) {
		double[] probd = readProbabilities(variant, sampleIdx, genotype);
		double dosage = probd[1] + 2 * probd[2];
		gts.add(new GenotypeBuilder(genotype).attribute("DS", outFormat.format(dosage)).make());	
	}

	private double[] readProbabilities(VariantContext variant, int sampleIdx, Genotype g) {
		String gp = (String)g.getAnyAttribute("GP");
		
		String[] probs = gp.split(",");
		if (probs.length != 3) {
			log.error("Illegal probabilities string for variant: " + variant.getID() + "; sample: " + sampleIdx);
			System.exit(97);
		}
		
		double[] probd = new double[3];
		try {
			probd[0] = Double.parseDouble(probs[0]);
			probd[1] = Double.parseDouble(probs[1]);
			probd[2] = Double.parseDouble(probs[2]);
		} catch (NumberFormatException nfe) {
			log.error("Not a number for genotype probability. Variant: " + variant.getID() + "; sample: " + sampleIdx);
			System.exit(96);
		}
		
		double sum = probd[0] + probd[1] + probd[2]; 
		if (Math.abs(sum - 1d) > 0.0021) { 
			// 3 significant digits, and give some tolerance for rounding
			log.error("Probabilities must sum up to 1 for variant: " + variant.getID() + "; sample: " + sampleIdx + "; sum = " + sum);
			System.exit(96);
		}

		return probd;
	}

	private void prepareWriteVcf() {
		// add header line for new "DS" field
		header.addMetaDataLine(new VCFFormatHeaderLine("DS", 1, VCFHeaderLineType.Float, "Genotype dosage"));
		
		// don't index as this needs dictionary
		VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
			.unsetOption(Options.INDEX_ON_THE_FLY); 
		writer = builder.setOutputFile(outputFile).build();
		
		writer.writeHeader(header);		
	}

	private void prepareReadVcf() {
		AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(inputFile.getAbsolutePath(), new VCFCodec(), false);

        header = (VCFHeader)reader.getHeader();
		if (!header.hasFormatLine("GP")) {
			System.err.println("File does not seem to have genotype probabilities.");
			System.exit(98);
		}
		
        progressLogger = new ProgressLogger(log, 10000);
	}

	private void parseCommandLine(String[] args) {
		if (args.length != 2) {
			log.error("Usage: AddVcfDosage <invcf> <outvcf> [<GP_field_name>]");
			System.exit(99);
		}
		
		inputFile = new File(args[0]);
		if (!inputFile.exists()) {
			log.error("Input file does not exist: " + args[0]);
			System.exit(99);
		}
		
		outputFile = new File(args[1]);
		if (outputFile.exists()) {
			log.error("Output file exists: " + args[1]);
			System.exit(99);
		}
	}
	
	private File inputFile;
	private File outputFile;
	
	private AbstractFeatureReader<VariantContext, LineIterator> reader;
	private VCFHeader header;
	private VariantContextWriter writer;

	private NumberFormat outFormat = new DecimalFormat("0.###", DecimalFormatSymbols.getInstance(Locale.US));
	private Log log = Log.getInstance(AddVcfDosage.class);

	private ProgressLogger progressLogger;
	
}
