package com.hartwig.hmftools.breakpointinspector;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class BreakPointInspectorApplication {

    private static final Logger LOGGER = LogManager.getLogger(BreakPointInspectorApplication.class);

    private static final String REF_BAM_PATH = "ref_bam";
    private static final String REF_BAM_SLICED_OUTPUT_PATH = "ref_bam_sliced_output";
    private static final String TUMOR_BAM_PATH = "tumor_bam";
    private static final String TUMOR_BAM_SLICED_OUTPUT_PATH = "tumor_bam_sliced_output";
    private static final String PROXIMITY = "proximity";
    private static final String VCF_IN = "input_vcf";
    private static final String VCF_OUT = "output_vcf";
    private static final String TSV_OUT = "output_tsv";
    private static final String CONTAMINATION_FRACTION = "contamination_fraction";

    public static void main(final String... args) throws IOException, ParseException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String refBamPath = cmd.getOptionValue(REF_BAM_PATH);
        final String refBamSlicedOutputPath = cmd.getOptionValue(REF_BAM_SLICED_OUTPUT_PATH);
        final String tumorBamPath = cmd.getOptionValue(TUMOR_BAM_PATH);
        final String tumorBamSlicedOutputPath = cmd.getOptionValue(TUMOR_BAM_SLICED_OUTPUT_PATH);
        final String vcfInputPath = cmd.getOptionValue(VCF_IN);
        final String vcfOutputPath = cmd.getOptionValue(VCF_OUT);
        final String tsvOutputPath = cmd.getOptionValue(TSV_OUT);

        if (refBamPath == null || tumorBamPath == null || vcfInputPath == null || tsvOutputPath == null) {
            printHelpAndExit(options);
            return;
        }

        final VCFFileReader vcfReader = new VCFFileReader(new File(vcfInputPath), false);

        final List<String> samples = vcfReader.getFileHeader().getGenotypeSamples();
        if (samples.size() != 2) {
            LOGGER.warn("Could not determine tumor and sample from VCF");
            System.exit(1);
        }

        final SamReader refReader = SamReaderFactory.makeDefault().open(new File(refBamPath));
        final SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(tumorBamPath));
        final Analysis analysis = buildAnalysis(cmd, refReader, tumorReader);

        LOGGER.info(String.format("Starting BPI filtering on vcf %s using tumor bam %s", vcfInputPath, tumorBamPath));
        final BPIAlgoOutput algo = BPIAlgo.run(vcfReader, tumorReader.getFileHeader().getSequenceDictionary(), analysis);
        LOGGER.info(String.format("Finishing BPI filtering. Generated %s variants", algo.variants().size()));

        if (vcfOutputPath != null) {
            writeToVCF(vcfOutputPath, vcfReader, algo.variants());
        }

        if (refBamSlicedOutputPath != null) {
            writeToSlice(refBamSlicedOutputPath, refReader, algo.optimizedIntervals());
        }

        if (tumorBamSlicedOutputPath != null) {
            writeToSlice(tumorBamSlicedOutputPath, tumorReader, algo.optimizedIntervals());
        }

        final List<String> tsv = Lists.newArrayList();
        tsv.add(TSVOutput.generateHeaders());
        tsv.addAll(algo.tsvOutput());
        Files.write(new File(tsvOutputPath).toPath(), tsv);

        refReader.close();
        tumorReader.close();
    }

    @NotNull
    private static Analysis buildAnalysis(@NotNull final CommandLine cmd, @NotNull final SamReader refReader,
            @NotNull final SamReader tumorReader) {
        final AnalysisBuilder analysisBuilder = new AnalysisBuilder();

        if (cmd.hasOption(PROXIMITY)) {
            analysisBuilder.range(Integer.parseInt(cmd.getOptionValue(PROXIMITY, "500")));
        }

        if (cmd.hasOption(CONTAMINATION_FRACTION)) {
            analysisBuilder.contaminationFraction(Float.parseFloat(cmd.getOptionValue(CONTAMINATION_FRACTION, "0")));
        }

        return analysisBuilder.refReader(refReader).tumorReader(tumorReader).createAnalysis();
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(Option.builder(REF_BAM_PATH).required().hasArg().desc("The reference BAM (required)").build());
        options.addOption(Option.builder(REF_BAM_SLICED_OUTPUT_PATH)
                .hasArg()
                .desc("The sliced reference BAM to output (optional)")
                .build());
        options.addOption(Option.builder(TUMOR_BAM_PATH).required().hasArg().desc("The tumor BAM (required)").build());
        options.addOption(Option.builder(TUMOR_BAM_SLICED_OUTPUT_PATH).hasArg().desc("The sliced Tumor BAM to output (optional)").build());
        options.addOption(Option.builder(PROXIMITY).hasArg().desc("Distance to scan around breakpoint (optional, default=500)").build());
        options.addOption(Option.builder(VCF_IN).required().hasArg().desc("Manta VCF file to inspect (required)").build());
        options.addOption(Option.builder(VCF_OUT).hasArg().desc("VCF output file (optional)").build());
        options.addOption(Option.builder(TSV_OUT).hasArg().desc("TSV output file (required)").build());
        options.addOption(Option.builder(CONTAMINATION_FRACTION)
                .hasArg()
                .desc("fraction of allowable normal support per tumor support read")
                .build());
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static void printHelpAndExit(final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Break-Point-Inspector", "A second layer of filtering on top of Manta", options, "", true);
        System.exit(1);
    }

    private static void writeToVCF(@NotNull String path, @NotNull VCFFileReader vcfReader, @NotNull List<VariantContext> variants) {
        final VCFHeader header = vcfReader.getFileHeader();
        header.addMetaDataLine(new VCFInfoHeaderLine("BPI_START", 1, VCFHeaderLineType.Integer, "BPI adjusted breakend location"));
        header.addMetaDataLine(new VCFInfoHeaderLine("BPI_END", 1, VCFHeaderLineType.Integer, "BPI adjusted breakend location"));
        header.addMetaDataLine(new VCFInfoHeaderLine("BPI_AMBIGUOUS",
                0,
                VCFHeaderLineType.Flag,
                "BPI could not determine the breakpoints, inspect manually"));
        header.addMetaDataLine(new VCFHeaderLine("bpiVersion",
                BreakPointInspectorApplication.class.getPackage().getImplementationVersion()));
        Filter.updateVCFHeader(header);
        AlleleFrequency.updateVCFHeader(header);

        List<VariantContext> variantsToInclude = Lists.newArrayList(variants);
        variantsToInclude.sort(new VariantContextComparator(header.getSequenceDictionary()));

        final VariantContextWriter writer = new VariantContextWriterBuilder().setReferenceDictionary(header.getSequenceDictionary())
                .setOutputFile(path)
                .build();
        writer.writeHeader(header);
        variantsToInclude.forEach(writer::add);
        writer.close();
    }

    private static void writeToSlice(@NotNull String path, @NotNull SamReader reader, @NotNull QueryInterval[] intervals) {
        final File outputBAM = new File(path);
        final SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(), true, outputBAM);
        final SAMRecordIterator iterator = reader.queryOverlapping(intervals);
        while (iterator.hasNext()) {
            writer.addAlignment(iterator.next());
        }
        iterator.close();
        writer.close();
    }
}
