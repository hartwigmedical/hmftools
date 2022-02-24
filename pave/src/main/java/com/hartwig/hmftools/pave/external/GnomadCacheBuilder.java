package com.hartwig.hmftools.pave.external;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Collectors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class GnomadCacheBuilder
{
    private final String mInputVcf;
    private final String mOutputDir;
    private final String mOutputId;
    private final String mSpecificChromosome;
    private final double mFreqThreshold;

    private static final String GNOMAD_FILE = "gnomad_file";
    private static final String SPECIFIC_CHROMOSOME = "specific_chr";
    private static final String FREQ_THRESHOLD = "freq_threshold";

    public static final String GNOMAD_FILE_ID = "gnomad_variants";

    public GnomadCacheBuilder(final CommandLine cmd)
    {
        mOutputDir = parseOutputDir(cmd);
        mInputVcf = cmd.getOptionValue(GNOMAD_FILE);
        mOutputId = cmd.getOptionValue(OUTPUT_ID);
        mSpecificChromosome = cmd.getOptionValue(SPECIFIC_CHROMOSOME);
        mFreqThreshold = Double.parseDouble(cmd.getOptionValue(FREQ_THRESHOLD, "0"));
    }

    public static String formFileId(final String dir, final String chromosome, final String outputId)
    {
        String outputFile = dir + GNOMAD_FILE_ID;

        if(chromosome != null)
            outputFile += "_chr" + chromosome;

        if(outputId != null)
            outputFile += "_" + outputId;

        outputFile += ".csv";
        return outputFile;
    }

    public void run()
    {
        if(mInputVcf == null || !Files.exists(Paths.get(mInputVcf)))
        {
            PV_LOGGER.error("missing input file, exiting");
            System.exit(1);
        }

        PV_LOGGER.info("parsing Gnomad file({}) specificChr({}) frequencyThreshold({})",
                mInputVcf, mSpecificChromosome, mFreqThreshold);

        String outputFile = formFileId(mOutputDir, mSpecificChromosome, mOutputId);

        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(
                    mInputVcf, new VCFCodec(), false);

            BufferedWriter writer = createBufferedWriter(outputFile, false);

            if(mSpecificChromosome.isEmpty())
                writer.write("Chromosome,");

            writer.write("Position,Ref,Alt,Frequency");
            writer.newLine();

            int itemCount = 0;
            int filteredCount = 0;
            int belowFreqCount = 0;

            for(VariantContext context : reader.iterator())
            {
                ++itemCount;

                if(itemCount > 0 && (itemCount % 100000) == 0)
                {
                    PV_LOGGER.debug("processed {} variants, filtered({}) belowFreq({}) current location({}:{})",
                            itemCount, filteredCount, belowFreqCount, context.getContig(), context.getStart());
                }

                if(context.isFiltered())
                {
                    ++filteredCount;
                    continue;
                }

                double frequency = context.getAttributeAsDouble("AF", 0);

                if(mFreqThreshold > 0 && frequency < mFreqThreshold)
                {
                    ++belowFreqCount;
                    continue;
                }

                String chromosome = context.getContig();
                int position = context.getStart();

                String ref = context.getReference().getBaseString();
                String alt = context.getAlternateAlleles().stream().map(Allele::toString).collect(Collectors.joining(","));

                if(mSpecificChromosome.isEmpty())
                    writer.write(String.format("%s,", chromosome));

                writer.write(String.format("%d,%s,%s,%.5f", position, ref, alt, frequency));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            PV_LOGGER.error(" failed to read gnomad file): {}", e.toString());
        }

        PV_LOGGER.info("Gnomad file parse complete");
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        options.addOption(GNOMAD_FILE, true, "Gnomad VCF input file");
        options.addOption(FREQ_THRESHOLD, true, "Population frequency (AF) threshold to write VCF entry");
        options.addOption(SPECIFIC_CHROMOSOME, true, "Produce file per chromosome");
        addOutputOptions(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);
        setLogLevel(cmd);

        GnomadCacheBuilder gnomadCacheBuilder = new GnomadCacheBuilder(cmd);
        gnomadCacheBuilder.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
