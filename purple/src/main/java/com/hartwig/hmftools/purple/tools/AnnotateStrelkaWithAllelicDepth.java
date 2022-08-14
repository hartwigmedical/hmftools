package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.File;

import com.hartwig.hmftools.common.variant.strelka.StrelkaAllelicDepth;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
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
    private static final String VCF_IN = "in";
    private static final String VCF_OUT = "out";

    private final String mInputVCF;
    private final String mOutputVCF;
    private final VCFHeader mHeader;
    private final VCFFileReader mVcfReader;
    private final VariantContextWriter mVcfWriter;

    public static void main(final String... args)
    {
        final Options options = createOptions();
        try (final AnnotateStrelkaWithAllelicDepth application = new AnnotateStrelkaWithAllelicDepth(options, args))
        {
            application.addAllelicDepthField();
        }
        catch (ParseException e)
        {
            PPL_LOGGER.warn(e);
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

        mInputVCF = cmd.getOptionValue(VCF_IN);
        mOutputVCF = cmd.getOptionValue(VCF_OUT);

        mVcfReader = new VCFFileReader(new File(mInputVCF), false);
        mHeader = generateOutputHeader(mVcfReader.getFileHeader());
        mVcfWriter = new VariantContextWriterBuilder().setOutputFile(mOutputVCF)
                .setReferenceDictionary(mHeader.getSequenceDictionary())
                .setIndexCreator(new TabixIndexCreator(mHeader.getSequenceDictionary(), new TabixFormat()))
                .setOption(htsjdk.variant.variantcontext.writer.Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .build();
    }

    private void addAllelicDepthField()
    {
        PPL_LOGGER.info("Reading input VCF: {}", mInputVCF);
        mVcfWriter.writeHeader(mHeader);
        for(final VariantContext variant : mVcfReader)
        {
            mVcfWriter.add(StrelkaAllelicDepth.enrich(variant));
        }
        PPL_LOGGER.info("Writing output VCF: {}", mOutputVCF);
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
        mVcfReader.close();
        mVcfWriter.close();
        PPL_LOGGER.info("Complete");
    }
}
