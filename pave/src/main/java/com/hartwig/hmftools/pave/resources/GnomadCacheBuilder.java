package com.hartwig.hmftools.pave.resources;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.variant.pon.GnomadCommon.GNOMAD_FILE_ID;
import static com.hartwig.hmftools.common.variant.pon.GnomadCommon.formFileId;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveConstants.APP_NAME;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

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


    public GnomadCacheBuilder(final ConfigBuilder configBuilder)
    {
        mOutputDir = parseOutputDir(configBuilder);
        mInputVcf = configBuilder.getValue(GNOMAD_FILE);
        mOutputId = configBuilder.getValue(OUTPUT_ID);
        mSpecificChromosome = configBuilder.getValue(SPECIFIC_CHROMOSOME, "");
        mFreqThreshold = configBuilder.getDecimal(FREQ_THRESHOLD);
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
            VcfFileReader reader = new VcfFileReader(mInputVcf);

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

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(GNOMAD_FILE, true, "Gnomad VCF input file");
        configBuilder.addDecimal(FREQ_THRESHOLD, "Population frequency (AF) threshold to write VCF entry", 0);
        configBuilder.addFlag(SPECIFIC_CHROMOSOME, "Produce file per chromosome");

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        GnomadCacheBuilder gnomadCacheBuilder = new GnomadCacheBuilder(configBuilder);
        gnomadCacheBuilder.run();
    }
}
