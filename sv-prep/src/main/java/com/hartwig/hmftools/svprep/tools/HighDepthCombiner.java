package com.hartwig.hmftools.svprep.tools;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.svprep.SvCommon.DELIM;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.tools.HighDepthConfig.HIGH_DEPTH_REGION_MAX_GAP;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.BlacklistLocations;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class HighDepthCombiner
{
    private final List<String> mInputFiles;
    private final Map<String,List<List<HighDepthRegion>>> mChrSampleHighDepthRegions;
    private final Map<String,List<CombinedRegion>> mChrCombinedRegions;
    private final BlacklistLocations mRefefenceBlacklist;
    private final int mMinSampleCount;
    private final BufferedWriter mWriter;

    private static final String INPUT_FILES = "input_files";
    private static final String OUTPUT_FILE = "output_file";
    private static final String REF_BLACKLIST_FILE = "ref_blacklist_file";
    private static final String MIN_SAMPLE_COUNT = "min_sample_count";

    private static final int DEFAULT_MIN_SAMPLE_COUNT = 4;
    private static final int MIN_REGION_LENGTH = 11;

    public HighDepthCombiner(final CommandLine cmd)
    {
        mInputFiles = Arrays.stream(cmd.getOptionValue(INPUT_FILES).split(",")).collect(Collectors.toList());
        mWriter = initialiseWriter(cmd.getOptionValue(OUTPUT_FILE));
        mMinSampleCount = Integer.parseInt(cmd.getOptionValue(MIN_SAMPLE_COUNT, String.valueOf(DEFAULT_MIN_SAMPLE_COUNT)));

        mChrSampleHighDepthRegions = Maps.newHashMap();
        mChrCombinedRegions = Maps.newHashMap();

        mRefefenceBlacklist = new BlacklistLocations(cmd.getOptionValue(REF_BLACKLIST_FILE));
    }

    public void run()
    {
        if(mInputFiles.isEmpty())
        {
            SV_LOGGER.error("no input files specified");
            System.exit(1);
        }

        SV_LOGGER.info("combing {} high depth region files", mChrSampleHighDepthRegions.size());

        loadSampleRegions();

        mergeRegions();

        writeCombinedResults();

        closeBufferedWriter(mWriter);

        SV_LOGGER.info("high depth region combination complete");
    }

    private void mergeRegions()
    {
        for(Map.Entry<String,List<List<HighDepthRegion>>> entry : mChrSampleHighDepthRegions.entrySet())
        {
            String chromosome = entry.getKey();
            List<CombinedRegion> combinedRegions = Lists.newArrayList();
            mChrCombinedRegions.put(chromosome, combinedRegions);

            SV_LOGGER.info("merging chromosome({})", chromosome);
            int sampleIndex = 1;

            for(List<HighDepthRegion> regions : entry.getValue())
            {
                SV_LOGGER.trace("merging sample({})", sampleIndex++);

                for(HighDepthRegion region : regions)
                {
                    int index = 0;
                    boolean matched = false;

                    while(index < combinedRegions.size())
                    {
                        CombinedRegion combinedRegion = combinedRegions.get(index);
                        if(positionsOverlap(region.Region.start(), region.Region.end(), combinedRegion.Start, combinedRegion.End))
                        {
                            matched = true;
                            combinedRegion.addBases(region.Region.start(), region.Region.end());
                            break;
                        }
                        else if(region.Region.end() < combinedRegion.Start)
                        {
                            break;
                        }

                        ++index;
                    }

                    if(!matched)
                    {
                        CombinedRegion combinedRegion = new CombinedRegion(region.Region.start(), region.Region.end());
                        combinedRegions.add(index, combinedRegion);
                    }
                }
            }
        }
    }

    private class CombinedRegion
    {
        public int Start;
        public int End;
        public List<Integer> Depth;

        public CombinedRegion(final int start, final int end)
        {
            Start = start;
            End = end;
            Depth = Lists.newArrayList();

            for(int pos = start; pos <= end; ++pos)
            {
                Depth.add(1);
            }
        }

        public int length() { return End - Start + 1; }

        public void addBases(int start, int end)
        {
            // first add overlapping bases
            for(int i = 0; i < length(); ++i)
            {
                int pos = Start + i;

                if(positionWithin(pos, start, end))
                    Depth.set(i, Depth.get(i) + 1);
            }

            // then add bases before
            for(int pos = start; pos < Start; ++pos)
            {
                Depth.add(0, 1);
            }

            // and bases afterwards
            for(int pos = End + 1; pos <= end; ++pos)
            {
                Depth.add(1);
            }

            Start = min(Start, start);
            End = max(End, end);
        }
    }

    private BufferedWriter initialiseWriter(final String filename)
    {
        SV_LOGGER.info("writing output to {}", filename);

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Chromosome\tPosStart\tPosEnd\tSamplesMin\tSamplesMax");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to initialise writer: {}", e.toString());
        }

        return null;
    }

    private List<HighDepthRegion> mergeChromosomeRegions(final String chromosome)
    {
        List<HighDepthRegion> highDepthRegions = Lists.newArrayList();

        List<CombinedRegion> combinedRegions = mChrCombinedRegions.get(chromosome.toString());
        if(combinedRegions == null)
            return highDepthRegions;

        for(CombinedRegion region : combinedRegions)
        {
            HighDepthRegion currentRegion = null;

            for(int i = 0; i < region.Depth.size(); ++i)
            {
                int position = region.Start + i;
                int baseDepth = region.Depth.get(i);

                if(baseDepth >= mMinSampleCount)
                {
                    if(currentRegion == null)
                    {
                        currentRegion = new HighDepthRegion(new ChrBaseRegion(chromosome.toString(), position, position));
                        currentRegion.DepthMin = baseDepth;
                        currentRegion.DepthMax = baseDepth;
                        highDepthRegions.add(currentRegion);
                    }
                    else
                    {
                        // extend the region
                        currentRegion.Region.setEnd(position);
                        currentRegion.DepthMax = max(currentRegion.DepthMax, baseDepth);
                    }
                }
                else
                {
                    if(currentRegion == null)
                        continue;

                    if(position - currentRegion.Region.end() < HIGH_DEPTH_REGION_MAX_GAP)
                        continue;

                    // end this region
                    currentRegion = null;
                }
            }
        }

        List<BaseRegion> referenceRegions = mRefefenceBlacklist.getRegions(chromosome);

        if(referenceRegions != null)
        {
            for(BaseRegion refRegion : referenceRegions)
            {
                // merge any adjacent regions
                int index = 0;
                boolean matched = false;
                while(index < highDepthRegions.size())
                {
                    HighDepthRegion region = highDepthRegions.get(index);

                    if(region.Region.start() > refRegion.end())
                        break;

                    if(positionsOverlap(region.Region.start(), region.Region.end(), refRegion.start(), refRegion.end()))
                    {
                        region.Region.setStart(min(region.Region.start(), refRegion.start()));
                        region.Region.setEnd(max(region.Region.end(), refRegion.end()));
                        matched = true;

                        // check if subsequent regions can now be merged in
                        int nextIndex = index + 1;
                        while(nextIndex < highDepthRegions.size())
                        {
                            HighDepthRegion nextRegion = highDepthRegions.get(nextIndex);

                            if(!positionsOverlap(nextRegion.Region.start(), nextRegion.Region.end(), refRegion.start(), refRegion.end()))
                                break;

                            highDepthRegions.remove(nextIndex);
                        }

                        break;
                    }
                    else
                    {
                        ++index;
                    }
                }

                if(!matched)
                    highDepthRegions.add(index, new HighDepthRegion(new ChrBaseRegion(chromosome, refRegion.start(), refRegion.end())));
            }
        }

        return highDepthRegions;
    }

    public void writeCombinedResults()
    {
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            List<HighDepthRegion> highDepthRegions = mergeChromosomeRegions(chromosome.toString());

            if(!highDepthRegions.isEmpty())
            {
                try
                {
                    for(HighDepthRegion region : highDepthRegions)
                    {
                        if(region.Region.baseLength() < MIN_REGION_LENGTH)
                            continue;

                        // write as a BED file, so note the -1 on the start
                        mWriter.write(format("%s\t%d\t%d\t%d\t%d",
                                region.Region.Chromosome, region.Region.start() - 1, region.Region.end(), region.DepthMin, region.DepthMax));
                        mWriter.newLine();
                    }
                }
                catch(IOException e)
                {
                    SV_LOGGER.error(" failed to write region: {}", e.toString());
                }
            }
        }
    }

    private void loadSampleRegions()
    {
        int totalRegions = 0;

        for(String filename : mInputFiles)
        {
            try
            {
                List<String> lines = Files.readAllLines(Paths.get(filename));
                lines.remove(0);

                Map<String,List<HighDepthRegion>> chrRegions = Maps.newHashMap();

                for(String line : lines)
                {
                    String[] values = line.split(DELIM, -1);

                    String chromosome = values[0];

                    List<HighDepthRegion> regions = chrRegions.get(chromosome);

                    if(regions == null)
                    {
                        regions = Lists.newArrayList();
                        chrRegions.put(chromosome, regions);
                    }

                    HighDepthRegion region = new HighDepthRegion(new ChrBaseRegion(
                            chromosome, Integer.parseInt(values[1]), Integer.parseInt(values[2])));
                    region.DepthMin = Integer.parseInt(values[3]);
                    region.DepthMax = values.length >= 5 ? Integer.parseInt(values[4]) : region.DepthMin;
                    regions.add(region);
                    ++totalRegions;
                }

                for(Map.Entry<String,List<HighDepthRegion>> entry : chrRegions.entrySet())
                {
                    List<List<HighDepthRegion>> sampleRegions = mChrSampleHighDepthRegions.get(entry.getKey());

                    if(sampleRegions == null)
                    {
                        sampleRegions = Lists.newArrayList();
                        mChrSampleHighDepthRegions.put(entry.getKey(), sampleRegions);
                    }

                    sampleRegions.add(entry.getValue());
                }
            }
            catch(IOException e)
            {
                SV_LOGGER.error("failed to read high-depth regions file: {}", e.toString());
            }
        }

        SV_LOGGER.info("loaded {} high-depth regions from {} files", totalRegions, mInputFiles.size());
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(INPUT_FILES, true, "Input files");
        options.addOption(OUTPUT_FILE, true, "Output file");
        options.addOption(REF_BLACKLIST_FILE, true, "Reference blacklist file to include");
        options.addOption(MIN_SAMPLE_COUNT, true, "Min sample count to produce region");
        addOutputOptions(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        HighDepthCombiner highDepthCombiner = new HighDepthCombiner(cmd);
        highDepthCombiner.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
