package com.hartwig.hmftools.purple.tools;

import java.io.File;

import com.hartwig.hmftools.common.variant.strelka.StrelkaAllelicDepth;

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
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class AnnotateStrelkaWithAllelicDepth implements AutoCloseable
{
    private static final Logger LOGGER = LogManager.getLogger(AnnotateStrelkaWithAllelicDepth.class);

    private static final String VCF_IN = "in";
    private static final String VCF_OUT = "out";

    private final String inputVCF;
    private final String outputVCF;
    private final VCFHeader header;
    private final VCFFileReader vcfReader;
    private final VariantContextWriter vcfWriter;

    public static void main(final String... args)
    {
        final Options options = createOptions();
        try (final AnnotateStrelkaWithAllelicDepth application = new AnnotateStrelkaWithAllelicDepth(options, args))
        {
            application.addAllelicDepthField();
        } catch (ParseException e)
        {
            LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("AnnotateStrelkaWithAllelicDepth", options);
            System.exit(1);
        }
    }

    private AnnotateStrelkaWithAllelicDepth(final Options options, final String... args) throws ParseException
    {
        final CommandLine cmd = createCommandLine(args, options);
        if(!cmd.hasOption(VCF_IN))
        {
            throw new ParseException(VCF_IN + " is a mandatory argument");
        }

        if(!cmd.hasOption(VCF_OUT))
        {
            throw new ParseException(VCF_OUT + " is a mandatory argument");
        }

        inputVCF = cmd.getOptionValue(VCF_IN);
        outputVCF = cmd.getOptionValue(VCF_OUT);

        vcfReader = new VCFFileReader(new File(inputVCF), false);
        header = generateOutputHeader(vcfReader.getFileHeader());
        vcfWriter = new VariantContextWriterBuilder().setOutputFile(outputVCF)
                .setReferenceDictionary(header.getSequenceDictionary())
                .setIndexCreator(new TabixIndexCreator(header.getSequenceDictionary(), new TabixFormat()))
                .setOption(htsjdk.variant.variantcontext.writer.Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .build();
    }

    private void addAllelicDepthField()
    {
        LOGGER.info("Reading input VCF: {}", inputVCF);
        vcfWriter.writeHeader(header);
        for(final VariantContext variant : vcfReader)
        {
            vcfWriter.add(StrelkaAllelicDepth.enrich(variant));
        }
        LOGGER.info("Writing output VCF: {}", outputVCF);
    }

    @NotNull
    private VCFHeader generateOutputHeader(@NotNull final VCFHeader template)
    {
        final VCFHeader outputVCFHeader = new VCFHeader(template.getMetaDataInInputOrder(), template.getGenotypeSamples());
        outputVCFHeader.addMetaDataLine(VCFStandardHeaderLines.getFormatLine("AD"));
        return outputVCFHeader;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(VCF_IN, true, "Input strelka VCF (.gz supported)");
        options.addOption(VCF_OUT, true, "Output VCF with AD field (.gz supported)");
        return options;
    }

    @Override
    public void close()
    {
        vcfReader.close();
        vcfWriter.close();
        LOGGER.info("Complete");
    }
}
