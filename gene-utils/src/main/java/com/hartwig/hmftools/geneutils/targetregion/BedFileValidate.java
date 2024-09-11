package com.hartwig.hmftools.geneutils.targetregion;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.bed.BedFileReader.loadBedFileChrMap;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.lowerChromosome;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.bed.BedFileReader;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class BedFileValidate
{
    private static final String INPUT_FILE = "input_file";
    private static final String OUTPUT_FILE = "output_file";

    public BedFileValidate(final ConfigBuilder configBuilder)
    {
        String inputFile = configBuilder.getValue(INPUT_FILE);
        String outputFile = configBuilder.getValue(OUTPUT_FILE);

        try
        {
            List<ChrBaseRegion> chrRegions = BedFileReader.loadBedFile(inputFile, false);

            GU_LOGGER.info("loaded {} regions from BED file({})", chrRegions.size(), inputFile);

            boolean fixFile = outputFile != null;

            validateRegions(chrRegions, fixFile);

            if(outputFile != null)
                writeBedRegions(chrRegions, outputFile);

            GU_LOGGER.info("BED validation complete");
        }
        catch(Exception e)
        {
            GU_LOGGER.error("failed to read BED file({}): {}", inputFile, e.toString());
            System.exit(1);
        }
    }

    public void validateRegions(final List<ChrBaseRegion> chrRegions, boolean fixFile)
    {
        // check ordering
        int orderingErrors = 0;
        int overlaps = 0;
        int invalidRegions = 0;

        for(int i = 0; i < chrRegions.size() - 1; ++i)
        {
            ChrBaseRegion region = chrRegions.get(i);

            if(!HumanChromosome.contains(region.chromosome()) || region.start() > region.end())
            {
                GU_LOGGER.debug("region({}) invalid", region);
                ++invalidRegions;
                continue;
            }

            ChrBaseRegion nextRegion = chrRegions.get(i + 1);

            if(region.chromosome().equals(nextRegion.chromosome()))
            {
                if(region.start() > nextRegion.start())
                {
                    GU_LOGGER.debug("region({}) ordered after next({})", region, nextRegion);
                    ++orderingErrors;
                }
                else if(region.end() >= nextRegion.start())
                {
                    GU_LOGGER.debug("region({}) overlaps next({})", region, nextRegion);
                    ++overlaps;
                }
            }
            else if(!lowerChromosome(region.chromosome(), nextRegion.chromosome()))
            {
                GU_LOGGER.debug("region({}) later chromosome than next({})", region, nextRegion);
                ++orderingErrors;
            }
        }

        GU_LOGGER.info("summary: invalid({}) mis-ordered({}) overlaps({})", invalidRegions, orderingErrors, overlaps);

        // order
        if(fixFile)
        {
            Collections.sort(chrRegions);
        }

        // merge regions, remove duplicates
        int mergedCount = 0;
        int removedCount = 0;

        int index = 0;
        while(index < chrRegions.size() - 1)
        {
            ChrBaseRegion region = chrRegions.get(index);

            if(!HumanChromosome.contains(region.chromosome()) || region.start() > region.end())
            {
                chrRegions.remove(index);
                ++removedCount;
                continue;
            }

            int nextIndex = index + 1;

            while(index < chrRegions.size())
            {
                ChrBaseRegion nextRegion = chrRegions.get(nextIndex);

                if(!region.chromosome().equals(nextRegion.chromosome()))
                    break;

                if(region.end() < nextRegion.start())
                    break;

                if(region.end() < nextRegion.end())
                    region.setEnd(nextRegion.end());

                chrRegions.remove(nextIndex);
                ++mergedCount;
            }

            ++index;
        }

        if(mergedCount > 0 || removedCount > 0)
        {
            GU_LOGGER.info("removed {} regions, merged {} regions", removedCount, mergedCount);
        }
    }

    private void writeBedRegions(final List<ChrBaseRegion> chrRegions, final String outputFile)
    {
        GU_LOGGER.info("writing {} regions to BED file({})", chrRegions.size(), outputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile);

            for(ChrBaseRegion region : chrRegions)
            {
                writer.write(format("%s\t%d\t%d",region.Chromosome, region.start() - 1, region.end()));
                writer.newLine();
            }

            writer.close();

        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write to BED file: {}", e.toString());
            System.exit(1);
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(INPUT_FILE, true, "Input BED file");
        configBuilder.addConfigItem(OUTPUT_FILE, false, "Output BED filename");
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        new BedFileValidate(configBuilder);
    }
}
