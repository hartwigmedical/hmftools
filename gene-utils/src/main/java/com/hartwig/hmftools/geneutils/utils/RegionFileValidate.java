package com.hartwig.hmftools.geneutils.utils;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.lowerChromosome;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class RegionFileValidate
{
    private static final String INPUT_FILE = "input_file";
    private static final String OUTPUT_FILE = "output_file";

    private String mOutputHeader;

    public RegionFileValidate(final ConfigBuilder configBuilder)
    {
        String inputFile = configBuilder.getValue(INPUT_FILE);
        String outputFile = configBuilder.getValue(OUTPUT_FILE);

        mOutputHeader = "";

        try
        {
            List<RegionData> regionDataList = loadRegionData(inputFile);

            GU_LOGGER.info("loaded {} regions from regions file({})", regionDataList.size(), inputFile);

            boolean fixFile = outputFile != null;

            validateRegions(regionDataList, fixFile);

            if(outputFile != null)
                writeBedRegions(regionDataList, outputFile);

            GU_LOGGER.info("BED validation complete");
        }
        catch(Exception e)
        {
            GU_LOGGER.error("failed to read BED file({}): {}", inputFile, e.toString());
            System.exit(1);
        }
    }

    private static boolean isBedFile(final String filename)
    {
        return filename.endsWith(".bed") || filename.endsWith(".bed.gz");
    }

    private List<RegionData> loadRegionData(final String filename)
    {
        List<RegionData> regionDataList = Lists.newArrayList();

        boolean isBedFile = isBedFile(filename);

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));

            int chrIndex = 0;
            int posStartIndex = 1;
            int posEndIndex = 2;

            // check for headers
            String header = lines.get(0);

            if(header.contains(FLD_CHROMOSOME) || header.contains(FLD_CHROMOSOME.toLowerCase()))
            {
                mOutputHeader = header;
                lines.remove(0);
            }

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                String chromosome = values[chrIndex];
                int posStart = Integer.parseInt(values[posStartIndex]);

                if(isBedFile)
                    ++posStart;

                int posEnd = Integer.parseInt(values[posEndIndex]);

                RegionData regionData = new RegionData(chromosome, posStart, posEnd);

                for(int i = 3; i < values.length; ++i)
                {
                    regionData.OtherData.add(values[i]);
                }

                regionDataList.add(regionData);
            }
        }
        catch(Exception e)
        {
            GU_LOGGER.error("failed to read file({}): {}", filename, e.toString());
            System.exit(1);
        }

        return regionDataList;
    }

    private class RegionData extends ChrBaseRegion
    {
        public final List<String> OtherData;

        public RegionData(final String chromosome, final int posStart, final int posEnd)
        {
            super(chromosome, posStart, posEnd);
            OtherData = Lists.newArrayList();
        }
    }

    public void validateRegions(final List<RegionData> regionDataList, boolean fixFile)
    {
        // check ordering
        int orderingErrors = 0;
        int overlaps = 0;
        int invalidRegions = 0;

        for(int i = 0; i < regionDataList.size() - 1; ++i)
        {
            ChrBaseRegion region = regionDataList.get(i);

            if(!HumanChromosome.contains(region.chromosome()) || region.start() > region.end())
            {
                GU_LOGGER.debug("region({}) invalid", region);
                ++invalidRegions;
                continue;
            }

            ChrBaseRegion nextRegion = regionDataList.get(i + 1);

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
            Collections.sort(regionDataList);
        }

        // merge regions, remove duplicates
        int mergedCount = 0;
        int removedCount = 0;

        int index = 0;
        while(index < regionDataList.size() - 1)
        {
            RegionData regionData = regionDataList.get(index);

            if(!HumanChromosome.contains(regionData.chromosome()) || regionData.start() > regionData.end())
            {
                regionDataList.remove(index);
                ++removedCount;
                continue;
            }

            int nextIndex = index + 1;

            while(index < regionDataList.size())
            {
                ChrBaseRegion nextRegion = regionDataList.get(nextIndex);

                if(!regionData.chromosome().equals(nextRegion.chromosome()))
                    break;

                if(regionData.end() < nextRegion.start())
                    break;

                // no attempt to merge other data fields
                if(regionData.end() < nextRegion.end())
                    regionData.setEnd(nextRegion.end());

                regionDataList.remove(nextIndex);
                ++mergedCount;
            }

            ++index;
        }

        if(mergedCount > 0 || removedCount > 0)
        {
            GU_LOGGER.info("removed {} regions, merged {} regions", removedCount, mergedCount);
        }
    }

    private void writeBedRegions(final List<RegionData> regionDataList, final String outputFile)
    {
        GU_LOGGER.info("writing {} regions to file({})", regionDataList.size(), outputFile);

        boolean isBedFile = isBedFile(outputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile);

            if(!mOutputHeader.isEmpty())
            {
                writer.write(mOutputHeader);
                writer.newLine();
            }

            for(RegionData regionData : regionDataList)
            {
                StringJoiner sj = new StringJoiner(TSV_DELIM);
                sj.add(regionData.Chromosome);

                int posStart = regionData.start();

                if(isBedFile)
                    --posStart;

                sj.add(String.valueOf(posStart));
                sj.add(String.valueOf(regionData.end()));

                regionData.OtherData.forEach(x -> sj.add(x));

                writer.write(sj.toString());
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

        configBuilder.addPath(INPUT_FILE, true, "Input regions TSV file - header is optional");
        configBuilder.addConfigItem(OUTPUT_FILE, false, "Output filename");
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        new RegionFileValidate(configBuilder);
    }
}
