package com.hartwig.hmftools.bamtools.depth;

import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.depth.HighDepthFinder.FLD_BASE_DEPTH_AVG;
import static com.hartwig.hmftools.bamtools.depth.HighDepthFinder.FLD_BASE_DEPTH_MAX;
import static com.hartwig.hmftools.bamtools.depth.HighDepthFinder.FLD_BASE_DEPTH_MIN;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.HighDepthRegion.FLD_DEPTH_AVG;
import static com.hartwig.hmftools.common.region.HighDepthRegion.FLD_DEPTH_MAX;
import static com.hartwig.hmftools.common.region.HighDepthRegion.FLD_DEPTH_MIN;
import static com.hartwig.hmftools.common.region.HighDepthRegion.FLD_SAMPLE_COUNT;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.HighDepthRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileDelimiters;

import org.jetbrains.annotations.NotNull;

public class HighDepthCombiner
{
    private final CombinerConfig mConfig;
    private final List<String> mInputFiles;

    private final GenicRegions mGenicRegions;

    private final Map<String,List<List<HighDepthRegion>>> mChrSampleHighDepthRegions; // per chromosome, per sample high-depth region data
    private final Map<String,List<HighDepthRegion>> mFinalRegions; // keyed by chromosome

    private static final int MIN_REGION_LENGTH = 11;
    protected static final int PANEL_HIGH_DEPTH_THRESHOLD = 2000;
    protected static final double CHROMOSOME_Y_SAMPLE_FRACTION = 0.4;

    public HighDepthCombiner(final ConfigBuilder configBuilder)
    {
        mConfig = new CombinerConfig(configBuilder);

        mInputFiles = Lists.newArrayList();

        for(String sampleId : mConfig.SampleIds)
        {
            mInputFiles.add(convertWildcardSamplePath(mConfig.HighDepthFiles, sampleId));
        }

        mGenicRegions = new GenicRegions(configBuilder);

        mChrSampleHighDepthRegions = Maps.newHashMap();
        mFinalRegions = Maps.newHashMap();
    }

    public void run()
    {
        if(mInputFiles.isEmpty())
        {
            BT_LOGGER.error("no input files specified");
            System.exit(1);
        }

        BT_LOGGER.info("combining {} high depth region files", mInputFiles.size());

        loadSampleRegions();

        List<CombinerMergeTask> mergeTasks = Lists.newArrayList();

        for(Map.Entry<String,List<List<HighDepthRegion>>> entry : mChrSampleHighDepthRegions.entrySet())
        {
            String chromosome = entry.getKey();
            List<List<HighDepthRegion>> sampleRegions = entry.getValue();

            mergeTasks.add(new CombinerMergeTask(mConfig, chromosome, sampleRegions));
        }

        List<Callable<Void>> callableList = mergeTasks.stream().collect(Collectors.toList());

        if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
            System.exit(1);

        for(CombinerMergeTask mergeTask : mergeTasks)
        {
            List<HighDepthRegion> highDepthRegions = mergeTask.highDepthRegions();

            if(!validateRegions(highDepthRegions))
                System.exit(1);

            mFinalRegions.put(mergeTask.chromosome(), highDepthRegions);
        }

        mGenicRegions.checkKnownGeneOverlaps(mFinalRegions);
        mGenicRegions.checkFixedGenicRegionOverlaps(mFinalRegions);

        writeCombinedResults();

        BT_LOGGER.info("High depth region combination complete");
    }

    private void writeCombinedResults()
    {
        BT_LOGGER.info("writing output to {}", mConfig.OutputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.OutputFile, false);

            StringJoiner header = new StringJoiner(TSV_DELIM);
            header.add(FLD_CHROMOSOME);
            header.add(FLD_POS_START);
            header.add(FLD_POS_END);

            if(mConfig.WriteWithLabel)
            {
                header.add("Label");
            }
            else
            {
                header.add(FLD_SAMPLE_COUNT);
                header.add(FLD_DEPTH_MIN);
                header.add(FLD_DEPTH_MAX);
                header.add(FLD_DEPTH_AVG);
            }

            writer.write(header.toString());
            writer.newLine();

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());
                List<HighDepthRegion> highDepthRegions = mFinalRegions.get(chrStr);

                if(highDepthRegions == null || highDepthRegions.isEmpty())
                    continue;

                for(HighDepthRegion region : highDepthRegions)
                {
                    if(region.baseLength() < MIN_REGION_LENGTH)
                        continue;

                    StringJoiner regionData = new StringJoiner(TSV_DELIM);
                    regionData.add(region.Chromosome);
                    regionData.add(String.valueOf(region.start() - 1));  // write as a BED file, so note the -1 on the start
                    regionData.add(String.valueOf(region.end()));

                    if(mConfig.WriteWithLabel)
                    {
                        regionData.add(String.format("HIGH_DEPTH_%d-%d_SC=%d", region.DepthMin, region.DepthMax, region.SampleCount));
                    }
                    else
                    {
                        regionData.add(String.valueOf(region.SampleCount));
                        regionData.add(String.valueOf(region.DepthMin));
                        regionData.add(String.valueOf(region.DepthMax));
                        regionData.add(String.valueOf(region.DepthAvg));
                    }

                    writer.write(regionData.toString());
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            BT_LOGGER.error(" failed to write final regions: {}", e.toString());
        }
    }

    private boolean validateRegions(final List<HighDepthRegion> regions)
    {
        for(int i = 0; i < regions.size() - 1; ++i)
        {
            HighDepthRegion region = regions.get(i);
            HighDepthRegion nextRegion = regions.get(i + 1);

            if(region.end() >= nextRegion.start())
            {
                BT_LOGGER.error("region({}) overlaps with next({})", region, nextRegion);
                return false;
            }
            else if(region.start() > nextRegion.start())
            {
                BT_LOGGER.error("region({}) after with next({})", region, nextRegion);
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
                String delim = FileDelimiters.inferFileDelimiter(filename);

                Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), delim);

                int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
                int posStartIndex = fieldsIndexMap.get(FLD_POS_START);
                int posEndIndex = fieldsIndexMap.get(FLD_POS_END);
                int depthMinIndex = fieldsIndexMap.get(FLD_BASE_DEPTH_MIN);
                int depthMaxIndex = fieldsIndexMap.get(FLD_BASE_DEPTH_MAX);
                Integer depthAvgIndex = fieldsIndexMap.get(FLD_BASE_DEPTH_AVG); // added in v1.7

                lines.remove(0);

                Map<String,List<HighDepthRegion>> chrRegions = Maps.newHashMap();

                for(String line : lines)
                {
                    String[] values = line.split(delim, -1);

                    String chromosome = values[chrIndex];

                    if(!mConfig.SpecificRegions.isEmpty() && mConfig.SpecificRegions.stream().noneMatch(x -> x.Chromosome.equals(chromosome)))
                        continue;

                    List<HighDepthRegion> regions = chrRegions.get(chromosome);

                    if(regions == null)
                    {
                        regions = Lists.newArrayList();
                        chrRegions.put(chromosome, regions);
                    }

                    int posStart = Integer.parseInt(values[posStartIndex]);
                    int posEnd = Integer.parseInt(values[posEndIndex]);

                    if(!mConfig.SpecificRegions.isEmpty() && mConfig.SpecificRegions.stream().noneMatch(x ->
                            x.Chromosome.equals(chromosome) && positionsOverlap(posStart, posEnd, x.start(), x.end())))
                        continue;

                    HighDepthRegion region = new HighDepthRegion(new ChrBaseRegion(chromosome, posStart, posEnd));
                    region.DepthMin = Integer.parseInt(values[depthMinIndex]);
                    region.DepthMax = Integer.parseInt(values[depthMaxIndex]);

                    if(depthAvgIndex != null)
                    {
                        region.DepthAvg = Integer.parseInt(values[depthAvgIndex]);

                        if(region.DepthAvg < region.DepthMin) // filter, now applied in the finder
                            continue;
                    }
                    else
                    {
                        // estimate from region length min and max - assumes steady increase from each edge then consistent depth at max
                        int readLength = max(mConfig.ReadLength, 100);
                        int regionLength = region.baseLength();
                        int estimatedMaxDepthLength = max(regionLength - 2 * readLength, 1);
                        double estimatedDepth = 0.5 * (regionLength + estimatedMaxDepthLength) * (region.DepthMax - region.DepthMin);
                        estimatedDepth += region.DepthMin * regionLength;
                        region.DepthAvg = (int)round(estimatedDepth / regionLength);
                    }

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
                BT_LOGGER.error("failed to read high-depth regions file: {}", e.toString());
            }
        }

        BT_LOGGER.info("loaded {} high-depth regions from {} files", totalRegions, mInputFiles.size());
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        CombinerConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        HighDepthCombiner highDepthCombiner = new HighDepthCombiner(configBuilder);
        highDepthCombiner.run();
    }
}
