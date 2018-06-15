package com.hartwig.hmftools.breakpointinspector;

import static java.util.Arrays.asList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiConsumer;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.breakpointinspector.datamodel.EnrichedVariantContext;
import com.hartwig.hmftools.breakpointinspector.datamodel.Range;

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
        final String tsvOutputPath = cmd.getOptionValue(TSV_OUT);

        if (refBamPath == null || tumorBamPath == null || vcfInputPath == null || tsvOutputPath == null) {
            printHelpAndExit(options);
            return;
        }

        final AnalysisBuilder analysisBuilder = new AnalysisBuilder();

        if (cmd.hasOption(PROXIMITY)) {
            analysisBuilder.range(Integer.parseInt(cmd.getOptionValue(PROXIMITY, "500")));
        }

        if (cmd.hasOption(CONTAMINATION_FRACTION)) {
            analysisBuilder.contaminationFraction(Float.parseFloat(cmd.getOptionValue(CONTAMINATION_FRACTION, "0")));
        }

        final SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(tumorBamPath));
        final SamReader refReader = SamReaderFactory.makeDefault().open(new File(refBamPath));
        final Analysis analysis = analysisBuilder.refReader(refReader).tumorReader(tumorReader).createAnalysis();

        final VCFFileReader vcfReader = new VCFFileReader(new File(vcfInputPath), false);

        final List<String> samples = vcfReader.getFileHeader().getGenotypeSamples();
        if (samples.size() != 2) {
            LOGGER.warn("Could not determine tumor and sample from VCF");
            System.exit(1);
        }

        final List<String> tsv = Lists.newArrayList();
        tsv.add(TSVOutput.generateHeaders());

        final List<QueryInterval> combinedQueryIntervals = Lists.newArrayList();

        final Map<String, VariantContext> variantPerIDMap = Maps.newHashMap();
        final List<VariantContext> variants = Lists.newArrayList();
        for (VariantContext variant : vcfReader) {
            variantPerIDMap.put(variant.getID(), variant);

            final VariantContext mateVariant = variant;
            if (variant.hasAttribute("MATEID")) {
                variant = variantPerIDMap.get(variant.getAttributeAsString("MATEID", ""));
                if (variant == null) {
                    continue;
                }
            }

            final EnrichedVariantContext enrichedVariant =
                    VariantEnrichment.enrich(variant, mateVariant, tumorReader.getFileHeader().getSequenceDictionary());

            final StructuralVariantResult result = analysis.processStructuralVariant(enrichedVariant);
            combinedQueryIntervals.addAll(asList(result.QueryIntervals));

            tsv.add(TSVOutput.generateVariant(variant, enrichedVariant, result));

            final BiConsumer<VariantContext, Boolean> vcfUpdater = (v, swap) -> {
                final Set<String> filters = v.getCommonInfo().getFiltersMaybeNull();
                if (filters != null) {
                    filters.clear();
                }
                // we will map BreakpointError to a flag
                if (result.Filters.contains(Filter.Filters.BreakpointError.toString())) {
                    v.getCommonInfo().putAttribute("BPI_AMBIGUOUS", true, true);
                } else {
                    v.getCommonInfo().addFilters(result.Filters);
                }
                if (result.Filters.isEmpty()) {
                    final List<Double> af = asList(result.AlleleFrequency.getLeft(), result.AlleleFrequency.getRight());
                    v.getCommonInfo().putAttribute(AlleleFrequency.VCF_INFO_TAG, swap ? Lists.reverse(af) : af, true);
                }
                if (result.Breakpoints.getLeft() != null) {
                    v.getCommonInfo().putAttribute(swap ? "BPI_END" : "BPI_START", result.Breakpoints.getLeft().Position, true);
                }
                if (result.Breakpoints.getRight() != null) {
                    v.getCommonInfo().putAttribute(swap ? "BPI_START" : "BPI_END", result.Breakpoints.getRight().Position, true);
                }

                // remove CIPOS / CIEND when we have an insert sequence
                if (!v.hasAttribute("IMPRECISE") && v.hasAttribute("SVINSSEQ")) {
                    v.getCommonInfo().removeAttribute("CIPOS");
                    v.getCommonInfo().removeAttribute("CIEND");
                }
                variants.add(v);
            };

            vcfUpdater.accept(variant, false);
            if (mateVariant != variant) {
                vcfUpdater.accept(mateVariant, true);
            }
        }

        // TODO: update START, END with BPI values and save Manta values in new attributes

        final String vcfOutputPath = cmd.getOptionValue(VCF_OUT);
        if (vcfOutputPath != null) {
            final VCFHeader header = vcfReader.getFileHeader();
            header.addMetaDataLine(new VCFInfoHeaderLine("BPI_START", 1, VCFHeaderLineType.Integer, "BPI adjusted breakend location"));
            header.addMetaDataLine(new VCFInfoHeaderLine("BPI_END", 1, VCFHeaderLineType.Integer, "BPI adjusted breakend location"));
            header.addMetaDataLine(new VCFInfoHeaderLine("BPI_AMBIGUOUS",
                    0,
                    VCFHeaderLineType.Flag,
                    "BPI could not determine the breakpoints, inspect manually"));
            header.addMetaDataLine(new VCFHeaderLine("bpiVersion",
                    BreakPointInspectorApplication.class.getPackage().getImplementationVersion()));
            Filter.UpdateVCFHeader(header);
            AlleleFrequency.UpdateVCFHeader(header);

            final VariantContextWriter writer = new VariantContextWriterBuilder().setReferenceDictionary(header.getSequenceDictionary())
                    .setOutputFile(vcfOutputPath)
                    .build();
            writer.writeHeader(header);
            variants.sort(new VariantContextComparator(header.getSequenceDictionary()));
            variants.forEach(writer::add);
            writer.close();
        }

        final QueryInterval[] optimizedIntervals =
                QueryInterval.optimizeIntervals(combinedQueryIntervals.toArray(new QueryInterval[combinedQueryIntervals.size()]));

        if (tumorBamSlicedOutputPath != null) {
            writeToSlice(tumorBamSlicedOutputPath, tumorReader, optimizedIntervals);
        }

        if (refBamSlicedOutputPath != null) {
            writeToSlice(refBamSlicedOutputPath, refReader, optimizedIntervals);
        }

        Files.write(new File(tsvOutputPath).toPath(), tsv);

        refReader.close();
        tumorReader.close();
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

    @NotNull
    private static Range fixup(@NotNull final Range uncertainty1, final boolean imprecise, final boolean inversion) {
        if (imprecise) {
            return uncertainty1;
        } else {
            return inversion ? Range.invert(uncertainty1) : uncertainty1;
        }
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
