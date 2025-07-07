package com.hartwig.hmftools.esvee.utils;

import static com.hartwig.hmftools.common.region.ChrBaseRegion.getChromosomeFieldIndex;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getPositionEndFieldIndex;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getPositionStartFieldIndex;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileDelimiters;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class BlacklistCompare
{
    private final Map<String,List<RegionInfo>> mOrigBlacklistLocations;
    private final Map<String,List<RegionInfo>> mNewBlacklistLocations;
    private final String mOutputFile;

    private static final String ORIG_BLACKLIST_FILE = "orig_blacklist_file";
    private static final String NEW_BLACKLIST_FILE = "new_blacklist_file";
    private static final String OUTPUT_FILE = "output_file";

    public BlacklistCompare(final ConfigBuilder configBuilder)
    {
        mOrigBlacklistLocations = loadBlackistFile(configBuilder.getValue(ORIG_BLACKLIST_FILE));
        mNewBlacklistLocations = loadBlackistFile(configBuilder.getValue(NEW_BLACKLIST_FILE));
        mOutputFile = configBuilder.getValue(OUTPUT_FILE);
    }

    public void run()
    {
        SV_LOGGER.info("Comparing blacklist files");

        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(FLD_CHROMOSOME);
            sj.add("OrigRegionStart").add("OrigRegionEnd");
            sj.add("NewRegionStart").add("NewRegionEnd");
            sj.add("OrigExtraInfo").add("NewExtraInfo");

            writer.write(sj.toString());
            writer.newLine();

            Set<String> chromosomes = mOrigBlacklistLocations.keySet().stream().collect(Collectors.toSet());
            mNewBlacklistLocations.keySet().forEach(x -> chromosomes.add(x));

            for(String chromosome : chromosomes)
            {
                List<RegionInfo> origRegions = mOrigBlacklistLocations.get(chromosome);
                List<RegionInfo> newRegions = mNewBlacklistLocations.get(chromosome);

                if(newRegions == null)
                {
                    for(RegionInfo regionInfo : origRegions)
                    {
                        writeComparison(chromosome, writer, regionInfo, null);
                    }
                }
                else if(origRegions == null)
                {
                    for(RegionInfo regionInfo : newRegions)
                    {
                        writeComparison(chromosome, writer, null, regionInfo);
                    }
                }
                else
                {
                    for(RegionInfo origRegion : origRegions)
                    {
                        RegionInfo newRegion = newRegions.stream().filter(x -> x.overlaps(origRegion)).findFirst().orElse(null);

                        if(newRegion != null)
                        {
                            writeComparison(chromosome, writer, origRegion, newRegion);
                            newRegions.remove(newRegion);
                        }
                        else
                        {
                            writeComparison(chromosome, writer, origRegion, null);
                        }
                    }

                    for(RegionInfo regionInfo : newRegions)
                    {
                        writeComparison(chromosome, writer, null, regionInfo);
                    }
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise output file: {}", e.toString());
            System.exit(1);
        }

        SV_LOGGER.info("Blacklist file comparison complete");
    }

    private void writeComparison(
            final String chromosome, final BufferedWriter writer,
            @Nullable final RegionInfo origRegion, @Nullable final RegionInfo newRegion) throws IOException
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(chromosome);

        if(origRegion != null)
        {
            sj.add(String.valueOf(origRegion.start()));
            sj.add(String.valueOf(origRegion.end()));
        }
        else
        {
            sj.add("N/A").add("N/A");
        }

        if(newRegion != null)
        {
            sj.add(String.valueOf(newRegion.start()));
            sj.add(String.valueOf(newRegion.end()));
        }
        else
        {
            sj.add("N/A").add("N/A");
        }

        sj.add(origRegion != null ? origRegion.ExtraInfo : "");
        sj.add(newRegion != null ? newRegion.ExtraInfo : "");

        writer.write(sj.toString());
        writer.newLine();
    }

    public Map<String,List<RegionInfo>> loadBlackistFile(final String filename)
    {
        Map<String,List<RegionInfo>> chrRegionsMap = Maps.newHashMap();

        try
        {
            String delim = FileDelimiters.inferFileDelimiter(filename);

            List<String> lines = Files.readAllLines(Paths.get(filename));
            boolean isBedFile = filename.contains(".bed");

            String header = lines.get(0);

            int chrIndex = 0;
            int posStartIndex = 1;
            int posEndIndex = 2;

            if(header.contains(FLD_CHROMOSOME))
            {
                Map<String, Integer> fieldIndexMap = FileReaderUtils.createFieldsIndexMap(lines.get(0), delim);
                lines.remove(0);

                chrIndex = getChromosomeFieldIndex(fieldIndexMap);
                posStartIndex = getPositionStartFieldIndex(fieldIndexMap);
                posEndIndex = getPositionEndFieldIndex(fieldIndexMap);
            }

            for(String line : lines)
            {
                final String[] values = line.split(delim, -1);

                String chromosome = values[chrIndex];

                List<RegionInfo> regions = chrRegionsMap.get(chromosome);

                if(regions == null)
                {
                    regions = Lists.newArrayList();
                    chrRegionsMap.put(chromosome, regions);
                }

                int posStart = Integer.parseInt(values[posStartIndex]);

                if(isBedFile)
                    ++posStart;

                int posEnd = Integer.parseInt(values[posEndIndex]);

                StringJoiner extraInfo = new StringJoiner(ITEM_DELIM);
                for(int i = posEndIndex + 1; i < values.length; ++i)
                {
                    extraInfo.add(values[i]);
                }

                regions.add(new RegionInfo(posStart, posEnd, extraInfo.toString()));
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to read blacklist file({}): {}", filename, e.toString());
            System.exit(1);
        }

        return chrRegionsMap;
    }

    private class RegionInfo extends BaseRegion
    {
        public final String ExtraInfo;

        public RegionInfo(final int posStart, final int posEnd, final String extraInfo)
        {
            super(posStart, posEnd);
            ExtraInfo = extraInfo;
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        addSampleIdFile(configBuilder, false);
        configBuilder.addPath(ORIG_BLACKLIST_FILE, true, "Original Blacklist file");
        configBuilder.addPath(NEW_BLACKLIST_FILE, true, "New Blacklist file");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output filename");

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        BlacklistCompare blacklistCompare = new BlacklistCompare(configBuilder);
        blacklistCompare.run();
    }
}
