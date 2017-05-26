package com.hartwig.hmftools.breakpointinspector;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import htsjdk.samtools.*;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import static com.hartwig.hmftools.breakpointinspector.Util.*;
import static com.hartwig.hmftools.breakpointinspector.Stats.*;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class BreakPointInspectorApplication {

    private static final String REF_PATH = "ref";
    private static final String TUMOR_PATH = "tumor";
    private static final String BREAK_POINT1 = "bp1";
    private static final String BREAK_POINT2 = "bp2";
    private static final String PROXIMITY = "proximity";
    private static final String SV_LEN = "svlen";
    private static final String VCF = "vcf";

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(REF_PATH, true, "the Reference BAM (indexed)");
        options.addOption(TUMOR_PATH, true, "the Tumor BAM (indexed)");
        options.addOption(BREAK_POINT1, true, "position of first break point in chrX:123456 format");
        options.addOption(BREAK_POINT2, true, "position of second break point in chrX:123456 format (optional)");
        options.addOption(PROXIMITY, true, "base distance around breakpoint");
        options.addOption(SV_LEN, true, "length of the SV to inspect (>0)");
        options.addOption(VCF, true, "VCF file to batch inspect (can be compressed)");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static void printHelpAndExit(final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Break-Point-Inspector", "Inspect structural variants", options, "", true);
        System.exit(1);
    }

    private static List<String> parseMantaPRSR(final Genotype genotype) {
        String pr = (String) genotype.getExtendedAttribute("PR", "0,0");
        String sr = (String) genotype.getExtendedAttribute("SR", "0,0");
        return Stream.concat(Arrays.stream(pr.split(",")), Arrays.stream(sr.split(","))).collect(Collectors.toList());
    }

    public static void main(final String... args) throws ParseException, IOException {

        final Options options = createOptions();
        try {
            final CommandLine cmd = createCommandLine(options, args);

            // grab arguments
            final String refBAM = cmd.getOptionValue(REF_PATH);
            final String tumorBAM = cmd.getOptionValue(TUMOR_PATH);
            final String bp1String = cmd.getOptionValue(BREAK_POINT1);
            final String bp2String = cmd.getOptionValue(BREAK_POINT2);
            final String vcfPath = cmd.getOptionValue(VCF);
            final int range = Integer.parseInt(cmd.getOptionValue(PROXIMITY, "500"));

            if (refBAM == null || tumorBAM == null || (vcfPath != null && bp1String != null) || (bp2String != null
                    && cmd.getOptionValue(SV_LEN) != null))
                printHelpAndExit(options);

            // load the files
            final File tumorFile = new File(tumorBAM);
            final SamReader tumorReader = SamReaderFactory.makeDefault().open(tumorFile);
            final File refFile = new File(refBAM);
            final SamReader refReader = SamReaderFactory.makeDefault().open(refFile);

            // output the header
            final ArrayList<String> header = Lists.newArrayList("ID", "MANTA_BP1", "MANTA_BP2", "MANTA_SVLEN",
                    "MANTA_REF_PR_NORMAL", "MANTA_REF_PR_SUPPORT", "MANTA_REF_SR_NORMAL", "MANTA_REF_SR_SUPPORT",
                    "MANTA_TUMOR_PR_NORMAL", "MANTA_TUMOR_PR_SUPPORT", "MANTA_TUMOR_SR_NORMAL",
                    "MANTA_TUMOR_SR_SUPPORT", "MANTA_HOMSEQ", "MANTA_INSSEQ");
            header.addAll(prefixList(Sample.GetHeader(), "REF_"));
            header.addAll(prefixList(Sample.GetHeader(), "TUMOR_"));
            header.add("TUMOR_CLIP_INFO");
            System.out.println(String.join("\t", header));

            if (vcfPath != null) {
                final File vcfFile = new File(vcfPath);
                final VCFFileReader vcfReader = new VCFFileReader(vcfFile, false);

                // work out the reference sample
                final List<String> samples = vcfReader.getFileHeader().getGenotypeSamples();
                final Predicate<String> isRef = s -> s.endsWith("R") || s.endsWith("BL");
                final String refSampleName = samples.stream().filter(isRef).findFirst().orElse(null);
                final String tumorSampleName = samples.stream().filter(
                        s -> s.endsWith("T") || !isRef.test(s)).findFirst().orElse(null);
                if (refSampleName == null || tumorSampleName == null) {
                    System.err.println("could not determine tumor and sample from VCF");
                    System.exit(1);
                    return;
                }

                for (final VariantContext variant : vcfReader) {
                    if (variant.isFiltered())
                        continue;

                    final String location = variant.getContig() + ":" + Integer.toString(variant.getStart());
                    final Location location1 = Location.parseLocationString(location,
                            tumorReader.getFileHeader().getSequenceDictionary());
                    final Location location2;

                    final String variantType = variant.getAttributeAsString("SVTYPE", "");
                    switch (variantType) {
                        case "DEL":
                        case "DUP":
                        case "INV":
                            final int svLen = Math.abs(variant.getAttributeAsInt("SVLEN", 0));
                            location2 = location1.add(svLen);
                            break;
                        case "BND":
                            final String call = variant.getAlternateAllele(0).getDisplayString();
                            final String[] split = call.split("[\\]\\[]");
                            location2 = Location.parseLocationString(split[1],
                                    tumorReader.getFileHeader().getSequenceDictionary());
                            break;
                        default:
                            System.err.println(variant.getID() + " : UNEXPECTED SVTYPE=" + variantType);
                            continue;
                    }

                    final List<String> extraData = Lists.newArrayList(variant.getID(), location1.toString(),
                            location2.toString(), variant.getAttributeAsString("SVLEN", "."));

                    extraData.addAll(parseMantaPRSR(variant.getGenotype(refSampleName)));
                    extraData.addAll(parseMantaPRSR(variant.getGenotype(tumorSampleName)));

                    extraData.add(variant.getAttributeAsString("HOMSEQ", "."));
                    extraData.add(variant.getAttributeAsString("SVINSSEQ", "."));

                    Analysis.processStructuralVariant(extraData, refReader, tumorReader, location1, location2, range);
                }
            } else {
                final int svLen = Integer.parseInt(cmd.getOptionValue(SV_LEN, "0"));

                // parse the location strings
                final Location location1 = Location.parseLocationString(bp1String,
                        tumorReader.getFileHeader().getSequenceDictionary());
                final Location location2;
                if (bp2String != null) {
                    location2 = Location.parseLocationString(bp2String,
                            tumorReader.getFileHeader().getSequenceDictionary());
                } else if (svLen > 0) {
                    location2 = location1.add(svLen);
                } else {
                    printHelpAndExit(options);
                    System.exit(1);
                    return;
                }

                final List<String> extraData = Lists.newArrayList("manual", location1.toString(),
                        location2.toString(), svLen > 0 ? Integer.toString(svLen) : ".");
                extraData.addAll(Collections.nCopies(10, "."));
                Analysis.processStructuralVariant(extraData, refReader, tumorReader, location1, location2, range);
            }

        } catch (ParseException e) {
            printHelpAndExit(options);
            System.exit(1);
        }
    }
}
