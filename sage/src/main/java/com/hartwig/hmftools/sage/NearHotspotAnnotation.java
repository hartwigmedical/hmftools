package com.hartwig.hmftools.sage;

import java.io.File;
import java.io.IOException;

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

import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class NearHotspotAnnotation {

    private static final Logger LOGGER = LogManager.getLogger(NearHotspotAnnotation.class);
    private static final String HOTSPOT_FLAG = "HOTSPOT";
    private static final String HOTSPOT_DESCRIPTION = "Site is at a known hotspot location";

    private static final String NEAR_HOTSPOT_FLAG = "NEAR_HOTSPOT";
    private static final String NEAR_HOTSPOT_DESCRIPTION = "Variant within " + HotspotEnrichment.DISTANCE + " bases of hotspot";

    private static final String OUT_PATH = "out";
    private static final String INPUT_VCF = "input_vcf";
    private static final String HOTSPOT = "hotspot";

    public static void main(String[] args) throws IOException, ParseException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String inputFilePath = cmd.getOptionValue(INPUT_VCF);
        final String outputFilePath = cmd.getOptionValue(OUT_PATH);
        final String hotspotPath = cmd.getOptionValue(HOTSPOT);

        if (outputFilePath == null || hotspotPath == null || inputFilePath == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("NearHotspotAnnotation", options);
            System.exit(1);
        }

        final HotspotEnrichment hotspotEnrichment = HotspotEnrichment.fromHotspotsFile(hotspotPath);
        final VCFFileReader reader = new VCFFileReader(new File(inputFilePath), false);
        final VCFHeader header = generateOutputHeader(reader.getFileHeader());
        final VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputFilePath)
                .setReferenceDictionary(header.getSequenceDictionary())
                .setIndexCreator(new TabixIndexCreator(header.getSequenceDictionary(), new TabixFormat()))
                .setOption(htsjdk.variant.variantcontext.writer.Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .build();
        writer.writeHeader(header);

        for (VariantContext context : reader) {
            if (!context.hasAttribute(HOTSPOT_FLAG)) {
                if (hotspotEnrichment.isOnHotspot(context)) {
                    context = new VariantContextBuilder(context).attribute(HOTSPOT_FLAG, true).make();
                } else if (hotspotEnrichment.isNearHotspot(context)) {
                    context = new VariantContextBuilder(context).attribute(NEAR_HOTSPOT_FLAG, true).make();
                }
            }
            writer.add(context);
        }

        reader.close();
        writer.close();

    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(OUT_PATH, true, "Output file.");
        options.addOption(INPUT_VCF, true, "Input file.");
        options.addOption(HOTSPOT, true, "Hotspot input file.");
        return options;
    }

    @NotNull
    private static VCFHeader generateOutputHeader(@NotNull final VCFHeader template) {
        final VCFHeader outputVCFHeader = new VCFHeader(template.getMetaDataInInputOrder(), template.getSampleNamesInOrder());
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(HOTSPOT_FLAG, 0, VCFHeaderLineType.Flag, HOTSPOT_DESCRIPTION));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(NEAR_HOTSPOT_FLAG, 0, VCFHeaderLineType.Flag, NEAR_HOTSPOT_DESCRIPTION));

        return outputVCFHeader;
    }

}
