package com.hartwig.hmftools.bamtools.depth;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.depth.GenicRegions.FIXED_GENE_REGIONS;
import static com.hartwig.hmftools.bamtools.depth.GenicRegions.REMOVE_GENE_OVERLAPS;
import static com.hartwig.hmftools.bamtools.depth.HighDepthConfig.HIGH_DEPTH_REGION_MAX_GAP;
import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.addGenePanelOption;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.addKnownFusionFileOption;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificRegions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
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
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.ExcludedRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.file.FileDelimiters;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class HighDepthCombiner
{
    private final List<String> mInputFiles;

    private final int mThreads;
    private final int mMinSampleCount;
    private final int mMinRegionSize;
    private final RefGenomeVersion mRefGenVersion;

    private final GenicRegions mGenicRegions;

    private final List<ChrBaseRegion> mSpecificRegions;

    private final Map<String,List<List<HighDepthRegion>>> mChrSampleHighDepthRegions;
    private final Map<String,List<HighDepthRegion>> mFinalRegions;
    private final String mOutputFile;
    private final boolean mWriteWithLabel;

    // config
    private static final String HIGH_DEPTH_FILES = "high_depth_files";
    private static final String OUTPUT_FILE = "output_file";
    private static final String MIN_SAMPLE_COUNT = "min_sample_count";
    private static final String MIN_REGION_SIZE = "min_region_size";
    private static final String WRITE_LABEL = "write_label";

    private static final int DEFAULT_MIN_SAMPLE_COUNT = 4;
    private static final int MIN_REGION_LENGTH = 11;
    protected static final int PANEL_HIGH_DEPTH_THRESHOLD = 2000;

    public HighDepthCombiner(final ConfigBuilder configBuilder)
    {
        List<String> sampleIds = loadSampleIdsFile(configBuilder);
        String highDepthFiles = configBuilder.getValue(HIGH_DEPTH_FILES);

        mInputFiles = Lists.newArrayList();

        for(String sampleId : sampleIds)
        {
            mInputFiles.add(convertWildcardSamplePath(highDepthFiles, sampleId));
        }

        mGenicRegions = new GenicRegions(configBuilder);

        mOutputFile = configBuilder.getValue(OUTPUT_FILE);
        mMinSampleCount = configBuilder.getInteger(MIN_SAMPLE_COUNT);
        mMinRegionSize = configBuilder.getInteger(MIN_REGION_SIZE);
        mThreads = parseThreads(configBuilder);
        mChrSampleHighDepthRegions = Maps.newHashMap();
        mFinalRegions = Maps.newHashMap();

        mRefGenVersion = RefGenomeVersion.from(configBuilder);

        mWriteWithLabel = configBuilder.hasFlag(WRITE_LABEL);

        mSpecificRegions = Lists.newArrayList();

        try
        {
            mSpecificRegions.addAll(loadSpecificRegions(configBuilder));
        }
        catch(ParseException e)
        {
            BT_LOGGER.error("failed to load specific regions");
        }
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

        List<MergeTask> mergeTasks = Lists.newArrayList();

        for(Map.Entry<String,List<List<HighDepthRegion>>> entry : mChrSampleHighDepthRegions.entrySet())
        {
            String chromosome = entry.getKey();
            List<List<HighDepthRegion>> sampleRegions = entry.getValue();

            mergeTasks.add(new MergeTask(chromosome, sampleRegions));
        }

        List<Callable<Void>> callableList = mergeTasks.stream().collect(Collectors.toList());

        if(!TaskExecutor.executeTasks(callableList, mThreads))
            System.exit(1);

        for(MergeTask mergeTask : mergeTasks)
        {
            List<HighDepthRegion> highDepthRegions = mergeTask.highDepthRegions();

            if(!validateRegions(highDepthRegions))
                System.exit(1);

            mFinalRegions.put(mergeTask.mChromosome, highDepthRegions);
        }

        mGenicRegions.checkKnownGeneOverlaps(mFinalRegions);
        mGenicRegions.checkFixedGenicRegionOverlaps(mFinalRegions);

        writeCombinedResults();

        BT_LOGGER.info("High depth region combination complete");
    }

    private class MergeTask implements Callable<Void>
    {
        private final String mChromosome;
        private List<List<HighDepthRegion>> mSampleRegions;
        private final List<CombinedRegion> mCombinedRegions;
        private final List<HighDepthRegion> mHighDepthRegions;

        public MergeTask(final String chromosome, final List<List<HighDepthRegion>> sampleRegions)
        {
            mChromosome = chromosome;
            mSampleRegions = sampleRegions;
            mCombinedRegions = Lists.newArrayList();
            mHighDepthRegions = Lists.newArrayList();
        }

        public String chromosome() { return mChromosome; }
        public List<HighDepthRegion> highDepthRegions() { return mHighDepthRegions; }

        public Void call()
        {
            // first generate regions from across all samples
            mergeSampleRegions();

            // then merge regions from all samples
            mHighDepthRegions.addAll(mergeChromosomeRegions(mCombinedRegions));
            return null;
        }

        private void mergeSampleRegions()
        {
            BT_LOGGER.info("merging chromosome({})", mChromosome);
            int sampleIndex = 1;

            for(List<HighDepthRegion> regions : mSampleRegions)
            {
                BT_LOGGER.trace("merging sample({})", sampleIndex++);

                for(HighDepthRegion region : regions)
                {
                    int index = 0;
                    boolean matched = false;

                    while(index < mCombinedRegions.size())
                    {
                        CombinedRegion combinedRegion = mCombinedRegions.get(index);
                        if(positionsOverlap(region.start(), region.end(), combinedRegion.start(), combinedRegion.end()))
                        {
                            matched = true;
                            combinedRegion.addBases(region);
                            break;
                        }
                        else if(region.end() < combinedRegion.start())
                        {
                            break;
                        }

                        ++index;
                    }

                    if(!matched)
                    {
                        CombinedRegion combinedRegion = new CombinedRegion(region);
                        mCombinedRegions.add(index, combinedRegion);
                    }
                    else
                    {
                        // check of this matched region now overlaps with following ones
                        CombinedRegion matchedRegion = mCombinedRegions.get(index);

                        int nextIndex = index + 1;
                        while(nextIndex < mCombinedRegions.size())
                        {
                            CombinedRegion combinedRegion = mCombinedRegions.get(nextIndex);

                            if(!positionsOverlap(matchedRegion.start(), matchedRegion.end(), combinedRegion.start(), combinedRegion.end()))
                                break;

                            matchedRegion.addRegion(combinedRegion);
                            mCombinedRegions.remove(nextIndex);
                        }
                    }
                }
            }
        }

        private List<HighDepthRegion> mergeChromosomeRegions(final List<CombinedRegion> combinedRegions)
        {
            List<HighDepthRegion> highDepthRegions = Lists.newArrayList();

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
                            currentRegion = new HighDepthRegion(new ChrBaseRegion(mChromosome, positionCount.Position, positionCount.Position));
                            currentRegion.DepthMin = positionCount.DepthMin;
                            currentRegion.DepthMax = positionCount.DepthMax;
                            currentRegion.SampleCount = positionCount.Count;
                            highDepthRegions.add(currentRegion);
                        }
                        else
                        {
                            // extend the region
                            currentRegion.setEnd(positionCount.Position);
                            currentRegion.DepthMin = min(currentRegion.DepthMin, positionCount.DepthMin);
                            currentRegion.DepthMax = max(currentRegion.DepthMax, positionCount.DepthMax);
                            currentRegion.SampleCount = max(currentRegion.SampleCount, positionCount.Count);
                        }
                    }
                    else
                    {
                        if(currentRegion == null)
                            continue;

                        if(positionCount.Position - currentRegion.end() < HIGH_DEPTH_REGION_MAX_GAP)
                            continue;

                        // end this region
                        currentRegion = null;
                    }
                }
            }

            // include the excluded region
            ChrBaseRegion excludedRegion = ExcludedRegions.getPolyGRegion(mRefGenVersion);
            List<BaseRegion> referenceRegions = Lists.newArrayList();

            if(excludedRegion.Chromosome.equals(mChromosome))
                referenceRegions.add(new BaseRegion(excludedRegion.start(), excludedRegion.end()));

            for(BaseRegion refRegion : referenceRegions)
            {
                // merge any adjacent regions
                int index = 0;
                boolean matched = false;
                while(index < highDepthRegions.size())
                {
                    HighDepthRegion region = highDepthRegions.get(index);

                    if(region.start() > refRegion.end())
                        break;

                    if(positionsOverlap(region.start(), region.end(), refRegion.start(), refRegion.end()))
                    {
                        matched = true;

                        // check if subsequent regions can now be merged in - and average out their min and max depth
                        long depthMinTotal = (long)region.baseLength() * region.DepthMin;
                        long depthMaxTotal = (long)region.baseLength() * region.DepthMax;
                        int regionBaseTotal = region.baseLength();

                        int nextIndex = index + 1;
                        while(nextIndex < highDepthRegions.size())
                        {
                            HighDepthRegion nextRegion = highDepthRegions.get(nextIndex);

                            if(!positionsOverlap(nextRegion.start(), nextRegion.end(), refRegion.start(), refRegion.end()))
                                break;

                            depthMinTotal += (long)nextRegion.baseLength() * nextRegion.DepthMin;
                            depthMaxTotal += (long)nextRegion.baseLength() * nextRegion.DepthMax;
                            regionBaseTotal += nextRegion.baseLength();


                            highDepthRegions.remove(nextIndex);
                        }

                        region.setStart(min(region.start(), refRegion.start()));
                        region.setEnd(max(region.end(), refRegion.end()));
                        region.DepthMin = (int)round(depthMinTotal / (double)regionBaseTotal);
                        region.DepthMax = (int)round(depthMaxTotal / (double)regionBaseTotal);
                        break;
                    }
                    else
                    {
                        ++index;
                    }
                }

                if(!matched)
                    highDepthRegions.add(index, new HighDepthRegion(new ChrBaseRegion(mChromosome, refRegion.start(), refRegion.end())));
            }

            // check min width for the region
            if(mMinRegionSize > 0)
            {
                for(HighDepthRegion highDepthRegion : highDepthRegions)
                {
                    if(highDepthRegion.baseLength() < mMinRegionSize)
                    {
                        int diff = mMinRegionSize - highDepthRegion.baseLength();
                        int halfExtension = diff / 2;

                        highDepthRegion.setStart(highDepthRegion.start() - halfExtension);
                        highDepthRegion.setEnd(highDepthRegion.end() + halfExtension);
                    }
                }

                ChrBaseRegion.checkMergeOverlaps(highDepthRegions, true);
            }

            return highDepthRegions;
        }
    }

    private void writeCombinedResults()
    {
        BT_LOGGER.info("writing output to {}", mOutputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            StringJoiner header = new StringJoiner(TSV_DELIM);
            header.add("Chromosome");
            header.add("PosStart");
            header.add("PosEnd");

            if(mWriteWithLabel)
            {
                header.add("Label");
            }
            else
            {
                header.add("SampleCount");
                header.add("DepthMin");
                header.add("DepthMax");
            }

            writer.write(header.toString());
            writer.newLine();

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrStr = mRefGenVersion.versionedChromosome(chromosome.toString());
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

                    if(mWriteWithLabel)
                    {
                        regionData.add(String.format("HIGH_DEPTH_%d-%d_SC=%d", region.DepthMin, region.DepthMax, region.SampleCount));
                    }
                    else
                    {
                        regionData.add(String.valueOf(region.SampleCount));
                        regionData.add(String.valueOf(region.DepthMin));
                        regionData.add(String.valueOf(region.DepthMax));
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

                lines.remove(0);

                Map<String,List<HighDepthRegion>> chrRegions = Maps.newHashMap();

                for(String line : lines)
                {
                    String[] values = line.split(delim, -1);

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
                BT_LOGGER.error("failed to read high-depth regions file: {}", e.toString());
            }
        }

        BT_LOGGER.info("loaded {} high-depth regions from {} files", totalRegions, mInputFiles.size());
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        addSampleIdFile(configBuilder, true);
        configBuilder.addConfigItem(HIGH_DEPTH_FILES, true, "High depth sample file(s), use '*' in for sampleId");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output file");
        configBuilder.addPath(FIXED_GENE_REGIONS, false, "Reference blacklist file to include");
        configBuilder.addInteger(MIN_SAMPLE_COUNT, "Min sample count to produce region", DEFAULT_MIN_SAMPLE_COUNT);
        configBuilder.addInteger(MIN_REGION_SIZE, "Min final region width", 0);
        configBuilder.addFlag(REMOVE_GENE_OVERLAPS, "Remove high depth regions that overlap driver or fusion genes");
        configBuilder.addFlag(WRITE_LABEL, "Write depth info as 'Label' column for compatibility with panel definition");
        configBuilder.addConfigItem(REF_GENOME_VERSION, REF_GENOME_VERSION_CFG_DESC);

        addGenePanelOption(configBuilder, false);
        addKnownFusionFileOption(configBuilder);
        addEnsemblDir(configBuilder);
        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        HighDepthCombiner highDepthCombiner = new HighDepthCombiner(configBuilder);
        highDepthCombiner.run();
    }
}
