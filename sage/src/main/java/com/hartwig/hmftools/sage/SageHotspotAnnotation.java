package com.hartwig.hmftools.sage;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.List;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SageHotspotAnnotation {

    private static final Logger LOGGER = LogManager.getLogger(SageHotspotApplication.class);

    // VCF FIELDS
    private static final String HOTSPOT_FLAG = "HOTSPOT";
    private static final String HOTSPOT_DESCRIPTION = "Site is at a known hotspot location";
    private static final String NEAR_HOTSPOT_FLAG = "NEAR_HOTSPOT";
    private static final String NEAR_HOTSPOT_DESCRIPTION = "Variant within " + HotspotEnrichment.DISTANCE + " bases of hotspot";
    private static final String RECOVERED_FLAG = "RECOVERED";
    private static final String RECOVERED_FLAG_DESCRIPTION = "Variant has been recovered";

    // Arguments
    private static final String OUT_VCF = "out";
    private static final String SOURCE_VCF = "source_vcf";
    private static final String HOTSPOT_VCF = "hotspot_vcf";
    private static final String KNOWN_HOTSPOTS = "known_hotspots";

    // Members
    private final HotspotEnrichment hotspotEnrichment;

    public static void main(String[] args) throws IOException, ParseException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String inputFilePath = cmd.getOptionValue(SOURCE_VCF);
        final String hotspotFilePath = cmd.getOptionValue(HOTSPOT_VCF);
        final String outputFilePath = cmd.getOptionValue(OUT_VCF);
        final String hotspotPath = cmd.getOptionValue(KNOWN_HOTSPOTS);

        if (outputFilePath == null || hotspotPath == null || inputFilePath == null || hotspotFilePath == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SageHotspotAnnotation", options);
            System.exit(1);
        }

        new SageHotspotAnnotation(hotspotPath).merge(inputFilePath, hotspotFilePath, outputFilePath);
    }

    private SageHotspotAnnotation(final String knownHotspotLocations) throws IOException {
        hotspotEnrichment = HotspotEnrichment.fromHotspotsFile(knownHotspotLocations);
    }

    private void merge(final String inputVcf, final String hotspotVcf, final String outputVCF) {
        final VCFFileReader inputReader = new VCFFileReader(new File(inputVcf), true);
        final VCFFileReader hotspotReader = new VCFFileReader(new File(hotspotVcf), false);
        final TreeSet<VariantContext> tree =
                new TreeSet<>(new VCComparator(inputReader.getFileHeader().getSequenceDictionary()));

        LOGGER.info("Loading somatic variants from {}", inputVcf);
        try (CloseableIterator<VariantContext> inputIterator = inputReader.iterator()) {
            while (inputIterator.hasNext()) {
                tree.add(annotate(inputIterator.next()));
            }
        }

        LOGGER.info("Loading hotspot variants from {}", hotspotVcf);
        final List<VariantContext> hotspots =
                hotspotReader.iterator().stream().map(SageHotspotAnnotation::recovered).collect(Collectors.toList());
        for (VariantContext hotspot : hotspots) {
            if (hotspot.isNotFiltered()) {
                final List<VariantContext> overlapping = inputReader.query(hotspot.getContig(), hotspot.getStart(), hotspot.getEnd()).toList();
                final VariantContext primary = primary(hotspot, overlapping);
                if (primary.equals(hotspot)) {
                    overlapping.forEach(tree::remove);
                    tree.add(hotspot);
                }
            }
        }

        LOGGER.info("Writing output to {}", outputVCF);
        final VCFHeader header = generateOutputHeader(inputReader.getFileHeader());
        final VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputVCF)
                .setReferenceDictionary(header.getSequenceDictionary())
                .setIndexCreator(new TabixIndexCreator(header.getSequenceDictionary(), new TabixFormat()))
                .setOption(htsjdk.variant.variantcontext.writer.Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .build();
        writer.writeHeader(header);

        for (VariantContext context : tree) {
            writer.add(context);
        }
        writer.close();
        ;

    }

    @NotNull
    private VariantContext annotate(@NotNull final VariantContext context) {
        VariantContextBuilder builder =
                new VariantContextBuilder(context).attribute(HOTSPOT_FLAG, false).attribute(NEAR_HOTSPOT_FLAG, false);

        if (hotspotEnrichment.isOnHotspot(context)) {
            return builder.attribute(HOTSPOT_FLAG, true).make();
        } else if (hotspotEnrichment.isNearHotspot(context)) {
            return builder.attribute(NEAR_HOTSPOT_FLAG, true).make();
        }

        return builder.make();
    }

    @NotNull
    private static VariantContext recovered(@NotNull final VariantContext context) {
        return new VariantContextBuilder(context).attribute(HOTSPOT_FLAG, true)
                .attribute(NEAR_HOTSPOT_FLAG, false)
                .attribute(RECOVERED_FLAG, true)
                .make();
    }

    @NotNull
    private static VariantContext primary(@NotNull final VariantContext hotspot, @NotNull final List<VariantContext> overlaps) {
        Comparator<VariantContext> comparator = Comparator.comparingInt(o -> o.getEnd() - o.getStart());
        final List<VariantContext> all = Lists.newArrayList(overlaps);
        all.add(hotspot);
        all.sort(comparator.reversed());

        return all.get(0);
    }

    @NotNull
    private static VCFHeader generateOutputHeader(@NotNull final VCFHeader template) {
        final VCFHeader outputVCFHeader = new VCFHeader(template.getMetaDataInInputOrder(), template.getSampleNamesInOrder());
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(HOTSPOT_FLAG, 0, VCFHeaderLineType.Flag, HOTSPOT_DESCRIPTION));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(NEAR_HOTSPOT_FLAG, 0, VCFHeaderLineType.Flag, NEAR_HOTSPOT_DESCRIPTION));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(RECOVERED_FLAG, 0, VCFHeaderLineType.Flag, RECOVERED_FLAG_DESCRIPTION));

        return outputVCFHeader;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(OUT_VCF, true, "Output file.");
        options.addOption(SOURCE_VCF, true, "Input file.");
        options.addOption(HOTSPOT_VCF, true, "Hotspot input file.");
        options.addOption(KNOWN_HOTSPOTS, true, "Tab separated file of known hotspot locations.");
        return options;
    }

    private class VCComparator extends VariantContextComparator {

        VCComparator(final SAMSequenceDictionary dictionary) {
            super(dictionary);
        }

        @Override
        public int compare(final VariantContext firstVariantContext, final VariantContext secondVariantContext) {
            int positionResult = super.compare(firstVariantContext, secondVariantContext);

            return positionResult == 0 ? firstVariantContext.getReference()
                    .getBaseString()
                    .compareTo(secondVariantContext.getReference().getBaseString()) : positionResult;
        }
    }

}
