package com.hartwig.hmftools.breakpointinspector;

import static java.util.Arrays.asList;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class BreakPointInspectorApplication {

    private static final String REF_PATH = "ref";
    private static final String REF_SLICE = "ref_slice";
    private static final String TUMOR_PATH = "tumor";
    private static final String TUMOR_SLICE = "tumor_slice";
    private static final String PROXIMITY = "proximity";
    private static final String VCF = "vcf";
    private static final String VCF_OUT = "output_vcf";

    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(Option.builder(REF_PATH).required().hasArg().desc("the Reference BAM (required)").build());
        options.addOption(Option.builder(REF_SLICE).hasArg().desc("the sliced Reference BAM to output (optional)").build());
        options.addOption(Option.builder(TUMOR_PATH).required().hasArg().desc("the Tumor BAM (required)").build());
        options.addOption(Option.builder(TUMOR_SLICE).hasArg().desc("the sliced Tumor BAM to output (optional)").build());
        options.addOption(Option.builder(PROXIMITY).hasArg().desc("distance to scan around breakpoint (optional, default=500)").build());
        options.addOption(Option.builder(VCF).required().hasArg().desc("Manta VCF file to batch inspect (required)").build());
        options.addOption(Option.builder(VCF_OUT).hasArg().desc("VCF output file (optional)").build());
        return options;
    }

    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static void printHelpAndExit(final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Break-Point-Inspector", "A second layer of filtering on top of Manta", options, "", true);
        System.exit(1);
    }

    public static void main(final String... args) throws ParseException, IOException {

        final AnalysisBuilder analysisBuilder = new AnalysisBuilder();
        final Options options = createOptions();
        try {
            final CommandLine cmd = createCommandLine(options, args);

            // grab arguments
            final String refPath = cmd.getOptionValue(REF_PATH);
            final String refSlicePath = cmd.getOptionValue(REF_SLICE);
            final String tumorPath = cmd.getOptionValue(TUMOR_PATH);
            final String tumorSlicePath = cmd.getOptionValue(TUMOR_SLICE);
            final String vcfPath = cmd.getOptionValue(VCF);

            if (cmd.hasOption(PROXIMITY)) {
                analysisBuilder.setRange(Integer.parseInt(cmd.getOptionValue(PROXIMITY, "500")));
            }

            if (refPath == null || tumorPath == null || vcfPath == null) {
                printHelpAndExit(options);
                return;
            }

            // load the files
            final File tumorBAM = new File(tumorPath);
            final SamReader tumorReader = SamReaderFactory.makeDefault().open(tumorBAM);
            final File refBAM = new File(refPath);
            final SamReader refReader = SamReaderFactory.makeDefault().open(refBAM);

            final File vcfFile = new File(vcfPath);
            final VCFFileReader vcfReader = new VCFFileReader(vcfFile, false);

            // get the sample names -- turns out Manta puts Ref first, then Tumor
            final List<String> samples = vcfReader.getFileHeader().getGenotypeSamples();
            if (samples.size() != 2) {
                System.err.println("could not determine tumor and sample from VCF");
                System.exit(1);
                return;
            }

            TSVOutput.PrintHeaders();

            final Analysis analysis = analysisBuilder.setRefReader(refReader).setTumorReader(tumorReader).createAnalysis();

            final List<QueryInterval> combinedQueryIntervals = Lists.newArrayList();

            final Map<String, VariantContext> variantMap = new HashMap<>();
            final List<VariantContext> variants = Lists.newArrayList();
            for (VariantContext variant : vcfReader) {

                variantMap.put(variant.getID(), variant);

                final VariantContext mateVariant = variant;
                if (variant.hasAttribute("MATEID")) {
                    variant = variantMap.get(variant.getAttributeAsString("MATEID", ""));
                    if (variant == null) {
                        continue;
                    }
                }

                final String location = variant.getContig() + ":" + Integer.toString(variant.getStart());
                Location location1 = Location.parseLocationString(location, tumorReader.getFileHeader().getSequenceDictionary());

                // uncertainty
                final List<Integer> CIPOS = variant.getAttributeAsIntList("CIPOS", 0);
                final Range uncertainty1 = CIPOS.size() == 2 ? new Range(CIPOS.get(0), CIPOS.get(1)) : new Range(0, 0);
                final List<Integer> CIEND = variant.getAttributeAsIntList("CIEND", 0);
                Range uncertainty2 = CIEND.size() == 2 ? new Range(CIEND.get(0), CIEND.get(1)) : new Range(0, 0);
                final boolean IMPRECISE = variant.hasAttribute("IMPRECISE");

                HMFVariantType svType;
                Location location2;
                switch (variant.getStructuralVariantType()) {
                    case INS:
                        svType = HMFVariantType.INS;
                        location2 = location1.set(variant.getAttributeAsInt("END", 0));
                        break;
                    case INV:
                        if (variant.hasAttribute("INV3")) {
                            svType = HMFVariantType.INV3;
                        } else if (variant.hasAttribute("INV5")) {
                            svType = HMFVariantType.INV5;
                        } else {
                            System.err.println(variant.getID() + " : expected either INV3 or INV5 flag");
                            continue;
                        }
                        location2 = location1.add(Math.abs(variant.getAttributeAsInt("SVLEN", 0)));
                        break;
                    case DEL:
                        svType = HMFVariantType.DEL;
                        location2 = location1.add(Math.abs(variant.getAttributeAsInt("SVLEN", 0)));
                        break;
                    case DUP:
                        svType = HMFVariantType.DUP;
                        location2 = location1.add(Math.abs(variant.getAttributeAsInt("SVLEN", 0)));
                        break;
                    case BND:
                        // process the breakend string
                        final String call = variant.getAlternateAllele(0).getDisplayString();
                        final String[] leftSplit = call.split("\\]");
                        final String[] rightSplit = call.split("\\[");

                        if (leftSplit.length >= 2) {
                            location2 = Location.parseLocationString(leftSplit[1], tumorReader.getFileHeader().getSequenceDictionary());
                            if (leftSplit[0].length() > 0) {
                                svType = HMFVariantType.INV3;
                                uncertainty2 = Range.invert(uncertainty1);
                            } else {
                                svType = HMFVariantType.DUP;
                                uncertainty2 = uncertainty1;
                            }
                        } else if (rightSplit.length >= 2) {
                            location2 = Location.parseLocationString(rightSplit[1], tumorReader.getFileHeader().getSequenceDictionary());
                            if (rightSplit[0].length() > 0) {
                                svType = HMFVariantType.DEL;
                                uncertainty2 = uncertainty1;
                            } else {
                                svType = HMFVariantType.INV5;
                                uncertainty2 = Range.invert(uncertainty1);
                            }
                        } else {
                            System.err.println(variant.getID() + " : could not parse breakpoint");
                            continue;
                        }

                        if (IMPRECISE) {
                            final List<Integer> MATE_CIPOS = mateVariant.getAttributeAsIntList("CIPOS", 0);
                            uncertainty2 = MATE_CIPOS.size() == 2 ? new Range(MATE_CIPOS.get(0), MATE_CIPOS.get(1)) : new Range(0, 0);
                        }

                        break;
                    default:
                        System.err.println(variant.getID() + " : UNEXPECTED SVTYPE=" + variant.getStructuralVariantType());
                        continue;
                }

                final HMFVariantContext ctx = new HMFVariantContext(variant.getID(), location1, location2, svType, IMPRECISE);
                ctx.Filter.addAll(variant.getFilters().stream().filter(s -> !s.startsWith("BPI")).collect(Collectors.toSet()));
                ctx.Uncertainty1 = uncertainty1;
                ctx.Uncertainty2 = uncertainty2;
                ctx.HomologySequence = variant.getAttributeAsString("HOMSEQ", "");
                ctx.InsertSequence = variant.getAttributeAsString("SVINSSEQ", "");
                ctx.BND = variant.getStructuralVariantType() == StructuralVariantType.BND;

                switch (ctx.Type) {
                    case INS:
                    case DEL:
                        ctx.OrientationBP1 = 1;
                        ctx.OrientationBP2 = -1;
                        break;
                    case INV3:
                        ctx.OrientationBP1 = 1;
                        ctx.OrientationBP2 = 1;
                        break;
                    case INV5:
                        ctx.OrientationBP1 = -1;
                        ctx.OrientationBP2 = -1;
                        break;
                    case DUP:
                        ctx.OrientationBP1 = -1;
                        ctx.OrientationBP2 = 1;
                        break;
                }

                final StructuralVariantResult result = analysis.processStructuralVariant(ctx);
                combinedQueryIntervals.addAll(asList(result.QueryIntervals));

                TSVOutput.PrintVariant(variant, ctx, result);

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
                // update header
                final VCFHeader header = vcfReader.getFileHeader();
                header.addMetaDataLine(new VCFInfoHeaderLine("BPI_START", 1, VCFHeaderLineType.Integer, "BPI adjusted breakend location"));
                header.addMetaDataLine(new VCFInfoHeaderLine("BPI_END", 1, VCFHeaderLineType.Integer, "BPI adjusted breakend location"));
                header.addMetaDataLine(new VCFInfoHeaderLine("BPI_AMBIGUOUS", 0, VCFHeaderLineType.Flag,
                        "BPI could not determine the breakends, inspect manually"));
                Filter.UpdateVCFHeader(header);
                AlleleFrequency.UpdateVCFHeader(header);
                // setup VCF
                final VariantContextWriter writer = new VariantContextWriterBuilder().setReferenceDictionary(header.getSequenceDictionary())
                        .setOutputFile(vcfOutputPath)
                        .build();
                writer.writeHeader(header);
                // write variants
                variants.sort(new VariantContextComparator(header.getSequenceDictionary()));
                variants.forEach(writer::add);
                // clean up
                writer.close();
            }

            // do a final slice pass

            final QueryInterval[] optimizedIntervals =
                    QueryInterval.optimizeIntervals(combinedQueryIntervals.toArray(new QueryInterval[combinedQueryIntervals.size()]));

            if (tumorSlicePath != null) {
                writeToSlice(tumorSlicePath, tumorReader, optimizedIntervals);
            }

            if (refSlicePath != null) {
                writeToSlice(refSlicePath, refReader, optimizedIntervals);
            }

            refReader.close();
            tumorReader.close();

        } catch (ParseException e) {
            printHelpAndExit(options);
            System.exit(1);
        }
    }

    private static void writeToSlice(final String path, final SamReader reader, final QueryInterval[] intervals) {
        final File tumorSliceBAM = new File(path);
        final SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(), true, tumorSliceBAM);
        final SAMRecordIterator iterator = reader.queryOverlapping(intervals);
        while (iterator.hasNext()) {
            writer.addAlignment(iterator.next());
        }
        iterator.close();
        writer.close();
    }
}
