package com.hartwig.hmftools.esvee.utils;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.isSingleBreakend;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.parseSingleOrientation;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.parseSvOrientation;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;
import static com.hartwig.hmftools.esvee.utils.BlacklistRepeatAnalyser.BASE_WINDOW_LENGTH;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.esvee.common.FilterType;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class RegionComplexityAnalyser
{
    private final List<String> mSampleIds;
    private final String mVcfFilename;
    private final String mOutputFile;

    private final int mBaseWindowLength;
    private final int mMinSamples;
    private final int mMinBreakends;

    private final int mThreads;
    private final SpecificRegions mSpecificChrRegions;

    private static final String VCF_FILE = "vcf_file";
    private static final String OUTPUT_FILE = "output_file";
    private static final String MIN_BREAKENDS = "min_breakends";
    private static final String MIN_SAMPLES = "min_samples";

    private static final int DEFAULT_BASE_WINDOW_LENGTH = 1000;

    public RegionComplexityAnalyser(final ConfigBuilder configBuilder)
    {
        mSampleIds = loadSampleIdsFile(configBuilder);
        mVcfFilename = configBuilder.getValue(VCF_FILE);
        mOutputFile = configBuilder.getValue(OUTPUT_FILE);

        mBaseWindowLength = configBuilder.getInteger(BASE_WINDOW_LENGTH);
        mMinBreakends = configBuilder.getInteger(MIN_BREAKENDS);
        mMinSamples = configBuilder.getInteger(MIN_SAMPLES);

        mThreads = parseThreads(configBuilder);
        mSpecificChrRegions = SpecificRegions.from(configBuilder);
    }

    public void run()
    {
        SV_LOGGER.info("analysing region complexity");

        long startTimeMs = System.currentTimeMillis();

        List<SampleTask> sampleTasks = Lists.newArrayList();
        for(int i = 0; i < min(mSampleIds.size(), mThreads); ++i)
        {
            sampleTasks.add(new SampleTask());
        }

        int taskIndex = 0;
        for(String sampleId : mSampleIds)
        {
            sampleTasks.get(taskIndex).sampleIds().add(sampleId);
            ++taskIndex;

            if(taskIndex >= sampleTasks.size())
                taskIndex = 0;
        }

        final List<Callable<Void>> callableList = sampleTasks.stream().collect(Collectors.toList());
        if(!TaskExecutor.executeTasks(callableList, mThreads))
            System.exit(1);

        // combined window data before writing
        WindowDataMap combinedWindowDataMap = new WindowDataMap();

        for(SampleTask sampleTask : sampleTasks)
        {
            combinedWindowDataMap.mergeOther(sampleTask.windowDataMap());
        }

        writeWindowData(combinedWindowDataMap);

        SV_LOGGER.info("Region complexity analysis complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private class SampleTask implements Callable<Void>
    {
        private final List<String> mSampleIds;
        private final WindowDataMap mWindowDataMap;

        public SampleTask()
        {
            mSampleIds = Lists.newArrayList();
            mWindowDataMap = new WindowDataMap();
        }

        public List<String> sampleIds()
        {
            return mSampleIds;
        }

        public WindowDataMap windowDataMap() { return mWindowDataMap; }

        @Override
        public Void call()
        {
            SV_LOGGER.info("processing {} samples", mSampleIds.size());

            int processed = 0;

            for(String sampleId : mSampleIds)
            {
                processSample(sampleId);

                ++processed;

                if((processed % 10) == 0)
                {
                    SV_LOGGER.debug("processed {} samples", processed);
                }
            }

            SV_LOGGER.debug("complete");

            return null;
        }

        private void processSample(final String sampleId)
        {
            String vcfFilename = mVcfFilename.replaceAll("\\*", sampleId);

            if(!Files.exists(Paths.get(vcfFilename)))
            {
                SV_LOGGER.warn("sample({}) vcf({}) doesn't exist, skipping", sampleId, vcfFilename);
                return;
            }

            VcfFileReader vcfFileReader = new VcfFileReader(vcfFilename);

            int processed = 0;
            String currentChromosome = "";
            WindowData currentWindow = null;
            int currentWindowBreakends = 0;

            for(VariantContext variantContext : vcfFileReader.iterator())
            {
                if(variantContext.getFilters().contains(FilterType.PON.vcfTag()))
                {
                    continue;
                }

                String chromosome = variantContext.getContig();
                int position = variantContext.getStart();

                if(mSpecificChrRegions.hasFilters() && mSpecificChrRegions.excludePosition(chromosome, position))
                    continue;

                byte orient = isSingleBreakend(variantContext) ? parseSingleOrientation(variantContext) : parseSvOrientation(variantContext);
                Orientation orientation = Orientation.fromByte(orient);

                int windowIndex = windowDataIndex(position);

                if(!chromosome.equals(currentChromosome) || currentWindow == null
                        || currentWindow.Index != windowIndex)
                {
                    currentChromosome = chromosome;

                    if(currentWindow != null)
                    {
                        currentWindow.addSampleBreakendCount(sampleId, currentWindowBreakends);
                    }

                    currentWindowBreakends = 0;
                    currentWindow = mWindowDataMap.getOrCreateWindow(chromosome, position);
                    ++currentWindow.SampleCount;
                }

                ++currentWindowBreakends;
                currentWindow.addBreakend(position, orientation);

                ++processed;

                if((processed % 100000) == 0)
                {
                    SV_LOGGER.debug("sample({}) processed {} variants", sampleId, processed);
                }
            }
        }
    }

    private class WindowData
    {
        public final int Index;
        public final Set<Integer> PosOrientationPositions;
        public final Set<Integer> NegOrientationPositions;

        public final List<SampleBreakendCount> TopSampleCounts;
        public int SampleCount;

        public WindowData(int index)
        {
            Index = index;
            PosOrientationPositions = Sets.newHashSet();
            NegOrientationPositions = Sets.newHashSet();

            TopSampleCounts = Lists.newArrayList();
            SampleCount = 0;
        }

        public void addBreakend(int position, final Orientation orientation)
        {
            if(orientation.isForward())
                PosOrientationPositions.add(position);
            else
                NegOrientationPositions.add(position);
        }

        private static final int MAX_SAMPLE_COUNT = 5;

        public void addSampleBreakendCount(final String sampleId, final int breakendCount)
        {
            if(TopSampleCounts.size() >= MAX_SAMPLE_COUNT && breakendCount <= TopSampleCounts.get(MAX_SAMPLE_COUNT - 1).BreakendCount)
                return;

            int index = 0;

            for(; index < TopSampleCounts.size(); ++index)
            {
                if(breakendCount > TopSampleCounts.get(index).BreakendCount)
                    break;
            }

            if(index >= MAX_SAMPLE_COUNT)
                return;

            TopSampleCounts.add(index, new SampleBreakendCount(sampleId, breakendCount));

            while(TopSampleCounts.size() > MAX_SAMPLE_COUNT)
            {
                TopSampleCounts.remove(TopSampleCounts.size() - 1);
            }
        }

        public void mergeOther(final WindowData other)
        {
            other.TopSampleCounts.forEach(x -> addSampleBreakendCount(x.SampleId, x.BreakendCount));

            other.PosOrientationPositions.forEach(x -> PosOrientationPositions.add(x));
            other.NegOrientationPositions.forEach(x -> NegOrientationPositions.add(x));
            SampleCount += other.SampleCount;
        }
    }

    private int windowDataIndex(final int position)
    {
        return position / mBaseWindowLength;
    }

    private class WindowDataMap
    {
        private final Map<String, Map<Integer, WindowData>> mChrWindowDataMap;

        public WindowDataMap()
        {
            mChrWindowDataMap = Maps.newHashMap();
        }

        public Map<String, Map<Integer, WindowData>> chrWindowDataMap() { return mChrWindowDataMap; }

        public WindowData getOrCreateWindow(final String chromosome, final int position)
        {
            Map<Integer, WindowData> windows = mChrWindowDataMap.get(chromosome);

            if(windows == null)
            {
                windows = Maps.newHashMap();
                mChrWindowDataMap.put(chromosome, windows);
            }

            int windowIndex = windowDataIndex(position);

            WindowData window = windows.get(windowIndex);

            if(window == null)
            {
                window = new WindowData(windowIndex);
                windows.put(windowIndex, window);
            }

            return window;
        }

        public void mergeOther(final WindowDataMap other)
        {
            for(Map.Entry<String, Map<Integer, WindowData>> chrEntry : other.chrWindowDataMap().entrySet())
            {
                String chromosome = chrEntry.getKey();

                for(Map.Entry<Integer, WindowData> windowEntry : chrEntry.getValue().entrySet())
                {
                    Map<Integer, WindowData> windows = mChrWindowDataMap.get(chromosome);

                    if(windows == null)
                    {
                        windows = Maps.newHashMap();
                        mChrWindowDataMap.put(chromosome, windows);
                    }

                    WindowData window = windows.get(windowEntry.getKey());

                    if(window == null)
                    {
                        window = new WindowData(windowEntry.getKey());
                        windows.put(windowEntry.getKey(), window);
                    }

                    window.mergeOther(windowEntry.getValue());
                }
            }
        }
    }

    private class SampleBreakendCount
    {
        public final String SampleId;
        public final int BreakendCount;

        public SampleBreakendCount(final String sampleId, final int breakendCount)
        {
            SampleId = sampleId;
            BreakendCount = breakendCount;
        }

        @Override public String toString() { return format("sampleId(%s) breakends(%d)", SampleId, BreakendCount); }
    }

    private void writeWindowData(final WindowDataMap windowDataMap)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(FLD_CHROMOSOME).add(FLD_POSITION_START).add(FLD_POSITION_END);

            sj.add("BreakendCount").add("SampleCount").add("MaxSamples");

            writer.write(sj.toString());
            writer.newLine();

            for(Map.Entry<String, Map<Integer, WindowData>> chrEntry : windowDataMap.chrWindowDataMap().entrySet())
            {
                String chromosome = chrEntry.getKey();

                for(WindowData window : chrEntry.getValue().values())
                {
                    if(window.SampleCount < mMinSamples)
                        continue;

                    int breakendCount = window.PosOrientationPositions.size() + window.NegOrientationPositions.size();

                    if(breakendCount < mMinBreakends)
                        continue;

                    sj = new StringJoiner(TSV_DELIM);

                    int windowPosStart = window.Index * mBaseWindowLength + 1;
                    int windowPosEnd = windowPosStart + mBaseWindowLength - 1;

                    sj.add(chromosome).add(String.valueOf(windowPosStart)).add(String.valueOf(windowPosEnd));
                    sj.add(String.valueOf(breakendCount));
                    sj.add(String.valueOf(window.SampleCount));

                    String sampleStr = window.TopSampleCounts.stream()
                            .map(x -> format("%s=%d", x.SampleId, x.BreakendCount)).collect(Collectors.joining(ITEM_DELIM));

                    sj.add(sampleStr);

                    writer.write(sj.toString());
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        addRefGenomeVersion(configBuilder);

        addSampleIdFile(configBuilder, false);
        configBuilder.addConfigItem(VCF_FILE, true, "VCF file, can use '*' in place of sampleIds");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output filename");

        configBuilder.addInteger(BASE_WINDOW_LENGTH, "Base window length for analysis", DEFAULT_BASE_WINDOW_LENGTH);
        configBuilder.addInteger(MIN_BREAKENDS, "Min breakends per window", 10);
        configBuilder.addInteger(MIN_SAMPLES, "Min samples per window", 3);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        RegionComplexityAnalyser regionComplexityAnalyser = new RegionComplexityAnalyser(configBuilder);
        regionComplexityAnalyser.run();
    }
}
