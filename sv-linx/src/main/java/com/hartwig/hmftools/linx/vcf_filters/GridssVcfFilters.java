package com.hartwig.hmftools.linx.vcf_filters;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class GridssVcfFilters
{
    private static final String GRIDSS_VCF_FILE = "vcf";
    private static final String SCOPE = "scope";
    private static final String SAMPLE = "sample";
    private static final String REF_GENOME = "ref_genome";
    private static final String OUTPUT_DIR = "output_dir";
    private static final String LOG_DEBUG = "log_debug";

    private static final Logger LOGGER = LogManager.getLogger(GridssVcfFilters.class);

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if(cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        final String sampleId = cmd.getOptionValue(SAMPLE);
        final String outputDir = cmd.getOptionValue(OUTPUT_DIR);

        final String vcfFile = cmd.getOptionValue(GRIDSS_VCF_FILE);
        final String refGenome = cmd.getOptionValue(REF_GENOME);
        final String scope = cmd.getOptionValue(SCOPE);

        LOGGER.info("scope({}) processing VCF({})", scope, vcfFile);

        loadVcf(vcfFile);

        LOGGER.info("VCF processing complete");
    }

    private static void loadVcf(final String vcfFile)
    {
        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                    vcfFile, new VCFCodec(), false);

            reader.iterator().forEach(x -> processVariant(x));
        }
        catch(Exception e)
        {
            LOGGER.error("error reading vcf({}): {}", vcfFile, e.toString());
        }
    }

    private static void processVariant(final VariantContext variant)
    {
        LOGGER.debug("position({}: {})", variant.getContig(), variant.getStart());

        // variant.

    }

    @NotNull
    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Name of the tumor sample");
        options.addOption(GRIDSS_VCF_FILE, true, "Path to the GRIDSS structural variant VCF file");
        options.addOption(REF_GENOME, true, "Ref genome fasta fila");
        options.addOption(SCOPE, true, "Scope: germilne or somatic");
        options.addOption(OUTPUT_DIR, true, "Path to write results");
        options.addOption(LOG_DEBUG, false, "Log verbose");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
