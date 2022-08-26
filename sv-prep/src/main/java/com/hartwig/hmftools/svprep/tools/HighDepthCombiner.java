package com.hartwig.hmftools.svprep.tools;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificRegions;
import static com.hartwig.hmftools.svprep.SvCommon.DELIM;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.tools.HighDepthConfig.HIGH_DEPTH_REGION_MAX_GAP;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.sv.ExcludedRegions;
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
    private final List<ChrBaseRegion> mSpecificRegions;
    private final RefGenomeVersion mRefGenVersion;
    private final BufferedWriter mWriter;

    private static final String HIGH_DEPTH_FILES = "high_depth_files";
    private static final String OUTPUT_FILE = "output_file";
    private static final String REF_BLACKLIST_FILE = "ref_blacklist_file";
    private static final String MIN_SAMPLE_COUNT = "min_sample_count";

    private static final int DEFAULT_MIN_SAMPLE_COUNT = 4;
    private static final int MIN_REGION_LENGTH = 11;

    public HighDepthCombiner(final CommandLine cmd)
    {
        List<String> sampleIds = loadSampleIdsFile(cmd);
        String highDepthFiles = cmd.getOptionValue(HIGH_DEPTH_FILES);

        mInputFiles = Lists.newArrayList();

        for(String sampleId : sampleIds)
        {
            mInputFiles.add(highDepthFiles.replaceAll("\\*", sampleId));
        }

        mWriter = initialiseWriter(cmd.getOptionValue(OUTPUT_FILE));
        mMinSampleCount = Integer.parseInt(cmd.getOptionValue(MIN_SAMPLE_COUNT, String.valueOf(DEFAULT_MIN_SAMPLE_COUNT)));

        mChrSampleHighDepthRegions = Maps.newHashMap();
        mChrCombinedRegions = Maps.newHashMap();

        mRefefenceBlacklist = new BlacklistLocations(cmd.getOptionValue(REF_BLACKLIST_FILE));
        mRefGenVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, V37.toString()));

        mSpecificRegions = Lists.newArrayList();

        try
        {
            mSpecificRegions.addAll(loadSpecificRegions(cmd));
        }
        catch(ParseException e)
        {
            SV_LOGGER.error("failed to load specific regions");
        }
    }

    public void run()
    {
        if(mInputFiles.isEmpty())
        {
            SV_LOGGER.error("no input files specified");
            System.exit(1);
        }

        SV_LOGGER.info("combining {} high depth region files", mInputFiles.size());

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
                        if(positionsOverlap(region.Region.start(), region.Region.end(), combinedRegion.start(), combinedRegion.end()))
                        {
                            matched = true;
                            combinedRegion.addBases(region.Region.start(), region.Region.end());
                            break;
                        }
                        else if(region.Region.end() < combinedRegion.start())
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
                    else
                    {
                        // check of this matched region now overlaps with following ones
                        CombinedRegion matchedRegion = combinedRegions.get(index);

                        int nextIndex = index + 1;
                        while(nextIndex < combinedRegions.size())
                        {
                            CombinedRegion combinedRegion = combinedRegions.get(nextIndex);

                            if(!positionsOverlap(matchedRegion.start(), matchedRegion.end(), combinedRegion.start(), combinedRegion.end()))
                                break;

                            matchedRegion.addRegion(combinedRegion);
                            combinedRegions.remove(nextIndex);
                        }
                    }
                }
            }
        }
    }

    private List<HighDepthRegion> mergeChromosomeRegions(final String chromosome)
    {
        List<HighDepthRegion> highDepthRegions = Lists.newArrayList();

        List<CombinedRegion> combinedRegions = mChrCombinedRegions.get(chromosome);
        if(combinedRegions == null)
            return highDepthRegions;

        for(CombinedRegion region : combinedRegions)
        {
            HighDepthRegion currentRegion = null;

            for(int i = 0; i < region.Depth.size(); ++i)
            {
                PositionCount positionCount = region.Depth.get(i);

                if(positionCount.Count >= mMinSampleCount)
                {
                    if(currentRegion == null)
                    {
                        currentRegion = new HighDepthRegion(new ChrBaseRegion(chromosome, positionCount.Position, positionCount.Position));
                        currentRegion.DepthMin = positionCount.Count;
                        currentRegion.DepthMax = positionCount.Count;
                        highDepthRegions.add(currentRegion);
                    }
                    else
                    {
                        // extend the region
                        currentRegion.Region.setEnd(positionCount.Position);
                        currentRegion.DepthMax = max(currentRegion.DepthMax, positionCount.Count);
                    }
                }
                else
                {
                    if(currentRegion == null)
                        continue;

                    if(positionCount.Position - currentRegion.Region.end() < HIGH_DEPTH_REGION_MAX_GAP)
                        continue;

                    // end this region
                    currentRegion = null;
                }
            }
        }

        List<BaseRegion> referenceRegions = mRefefenceBlacklist.getRegions(chromosome);

        // include the excluded region
        ChrBaseRegion excludedRegion = ExcludedRegions.getPolyGRegion(mRefGenVersion);

        if(excludedRegion.Chromosome.equals(chromosome))
        {
            referenceRegions.add(new BaseRegion(excludedRegion.start(), excludedRegion.end()));
        }

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

    private class CombinedRegion
    {
        public List<PositionCount> Depth;

        public CombinedRegion(final int start, final int end)
        {
            Depth = Lists.newArrayList();

            for(int pos = start; pos <= end; ++pos)
            {
                Depth.add(new PositionCount(pos));
            }
        }

        public int start() { return Depth.get(0).Position; }
        public int end() { return Depth.get(Depth.size() - 1).Position; }
        public int length() { return Depth.size(); }

        public void addBases(int start, int end)
        {
            for(int posNew = start; posNew <= end; ++posNew)
            {
                int existingIndex = 0;

                boolean found = false;
                while(existingIndex < Depth.size())
                {
                    PositionCount existing = Depth.get(existingIndex);

                    if(posNew == existing.Position)
                    {
                        found = true;
                        ++existing.Count;
                        break;
                    }

                    if(posNew < existing.Position)
                        break;

                    ++existingIndex;
                }

                if(!found)
                {
                    Depth.add(existingIndex, new PositionCount(posNew));
                }
            }
        }

        public void addRegion(final CombinedRegion other)
        {
            for(PositionCount otherCount : other.Depth)
            {
                int existingIndex = 0;

                boolean found = false;
                while(existingIndex < Depth.size())
                {
                    PositionCount existing = Depth.get(existingIndex);

                    if(otherCount.Position == existing.Position)
                    {
                        found = true;
                        existing.Count += otherCount.Count;
                        break;
                    }

                    if(otherCount.Position < existing.Position)
                        break;

                    ++existingIndex;
                }

                if(!found)
                {
                    Depth.add(existingIndex, otherCount);
                }
            }
        }

        public String toString() { return format("span(%d - %d) length(%d)", start(), end(), length()); }
    }

    private class PositionCount
    {
        public int Position;
        public int Count;

        public PositionCount(final int position)
        {
            Position = position;
            Count = 1;
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

    public void writeCombinedResults()
    {
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mRefGenVersion.versionedChromosome(chromosome.toString());
            List<HighDepthRegion> highDepthRegions = mergeChromosomeRegions(chrStr);

            if(!highDepthRegions.isEmpty())
            {
                if(!validateRegions(highDepthRegions))
                    System.exit(1);

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

    private boolean validateRegions(final List<HighDepthRegion> regions)
    {
        for(int i = 0; i < regions.size() - 1; ++i)
        {
            HighDepthRegion region = regions.get(i);
            HighDepthRegion nextRegion = regions.get(i + 1);

            if(region.Region.end() >= nextRegion.Region.start())
            {
                SV_LOGGER.error("region({}) overlaps with next({})", region, nextRegion);
                return false;
            }
            else if(region.Region.start() > nextRegion.Region.start())
            {
                SV_LOGGER.error("region({}) after with next({})", region, nextRegion);
                return false;
            }
        }

        return true;
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

                    if(!mSpecificRegions.isEmpty() && mSpecificRegions.stream().noneMatch(x -> x.Chromosome.equals(chromosome)))
                        continue;

                    List<HighDepthRegion> regions = chrRegions.get(chromosome);

                    if(regions == null)
                    {
                        regions = Lists.newArrayList();
                        chrRegions.put(chromosome, regions);
                    }

                    int posStart = Integer.parseInt(values[1]);
                    int posEnd = Integer.parseInt(values[2]);

                    if(!mSpecificRegions.isEmpty() && mSpecificRegions.stream().noneMatch(x ->
                            x.Chromosome.equals(chromosome) && positionsOverlap(posStart, posEnd, x.start(), x.end())))
                        continue;

                    HighDepthRegion region = new HighDepthRegion(new ChrBaseRegion(chromosome, posStart, posEnd));
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
        addSampleIdFile(options);
        options.addOption(HIGH_DEPTH_FILES, true, "High depth sample file(s), use '*' in for sampleId");
        options.addOption(OUTPUT_FILE, true, "Output file");
        options.addOption(REF_BLACKLIST_FILE, true, "Reference blacklist file to include");
        options.addOption(MIN_SAMPLE_COUNT, true, "Min sample count to produce region");
        options.addOption(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        addOutputOptions(options);
        addLoggingOptions(options);
        addSpecificChromosomesRegionsConfig(options);

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
