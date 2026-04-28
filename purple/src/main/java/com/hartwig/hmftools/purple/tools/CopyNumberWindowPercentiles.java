package com.hartwig.hmftools.purple.tools;

import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.CHR_PREFIX;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
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
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.purple.PurpleConstants;

import org.jetbrains.annotations.NotNull;

public class CopyNumberWindowPercentiles
{
    private final List<String> mSampleIds;
    private final String mPurpleDataDir;
    private final int mThreads;

    private final String mOutputFilename;
    private final boolean mCollapseOutput;

    private final Map<String,List<WindowValues>> mChrWindowsMap;

    private static final String OUTPUT_FILE = "output_file";
    private static final String COLLAPSE_OUTPUT = "collapse_output";

    private static final int WINDOW_SIZE = PurpleConstants.WINDOW_SIZE;
    private static final double TARGET_PERCENTILE_COPY_NUMBER = 1.0;

    public CopyNumberWindowPercentiles(final ConfigBuilder configBuilder)
    {
        mSampleIds = loadSampleIdsFile(configBuilder);
        mPurpleDataDir = configBuilder.getValue(PURPLE_DIR_CFG);
        mThreads = parseThreads(configBuilder);
        mOutputFilename = configBuilder.getValue(OUTPUT_FILE);
        mCollapseOutput = configBuilder.hasFlag(COLLAPSE_OUTPUT);

        mChrWindowsMap = Maps.newHashMap();
    }

    public void run()
    {
        if(mSampleIds.isEmpty())
        {
            PPL_LOGGER.error("missing sampleIds, exiting");
            System.exit(1);
        }

        PPL_LOGGER.info("generating cohort copy-number percentiles from {} samples", mSampleIds.size());
        long startTimeMs = System.currentTimeMillis();

        List<SampleCopyNumberTask> sampleTasks = Lists.newArrayList();

        for(int i = 0; i < min(mSampleIds.size(), mThreads); ++i)
        {
            sampleTasks.add(new SampleCopyNumberTask());
        }

        int taskIndex = 0;
        for(String sampleId : mSampleIds)
        {
            if(taskIndex >= sampleTasks.size())
                taskIndex = 0;

            sampleTasks.get(taskIndex).getSampleIds().add(sampleId);

            ++taskIndex;
        }

        List<Callable<Void>> callableList = sampleTasks.stream().collect(Collectors.toList());
        if(!TaskExecutor.executeTasks(callableList, mThreads))
        {
            System.exit(1);
        }

        PPL_LOGGER.info("writing cohort file: {}", mOutputFilename);

        writePercentilesFile();

        PPL_LOGGER.info("Purple copy-number percentiles generated, mins({})", runTimeMinsStr(startTimeMs));
    }

    private synchronized void addSampleCopyNumberValues(final String sampleId, final Map<String,List<Double>> chrCopyNumberValues)
    {
        for(Map.Entry<String,List<Double>> entry : chrCopyNumberValues.entrySet())
        {
            String chromosome = entry.getKey();
            List<Double> copyNumbers = entry.getValue();
            List<WindowValues> windows = mChrWindowsMap.get(chromosome);

            if(windows == null)
            {
                windows = Lists.newArrayListWithCapacity(copyNumbers.size());
                mChrWindowsMap.put(chromosome, windows);

                int positionStart = 1;
                for(double copyNumber : copyNumbers)
                {
                    WindowValues windowValues = new WindowValues(positionStart, mSampleIds.size());
                    windowValues.Values.add(copyNumber);
                    windows.add(windowValues);
                    positionStart += WINDOW_SIZE;
                }
            }
            else
            {
                for(int i = 0; i < copyNumbers.size(); ++i)
                {
                    if(windows.size() != copyNumbers.size())
                    {
                        PPL_LOGGER.error("sample({}) chromosome({}) mismatch windows({}) vs new CN({})",
                                sampleId, chromosome, windows.size(), copyNumbers.size());
                        System.exit(1);
                    }

                    WindowValues windowValues = windows.get(i);
                    windowValues.Values.add(copyNumbers.get(i));
                }
            }
        }
    }

    private class WindowValues
    {
        public final int PostitionStart;
        public final List<Double> Values;

        public WindowValues(final int postitionStart, final int sampleCount)
        {
            PostitionStart = postitionStart;
            Values = Lists.newArrayListWithCapacity(sampleCount);
        }

        public double percentile(final double targetCopyNumber)
        {
            if(Values.isEmpty())
                return -1;

            Collections.sort(Values);

            for(int i = 0; i < Values.size(); ++i)
            {
                double copyNumber = Values.get(i);
                if(abs(targetCopyNumber - copyNumber) < 0.001 || copyNumber > targetCopyNumber)
                {
                    return (i + 1) / (double)Values.size();
                }
            }

            return -1;
        }

        public void logValues(final String chromosome, int maxValues)
        {
            int valueCount = Values.size();

            // debug only
            if(Values.size() > maxValues * 2)
            {
                for(int i = 0; i < maxValues; ++i)
                {
                    PPL_LOGGER.debug(format("index(%d) location(%s:%d) copyNumber(%.3f)", i, chromosome, PostitionStart, Values.get(i)));
                }

                for(int i = valueCount - maxValues; i < valueCount; ++i)
                {
                    PPL_LOGGER.debug(format("index(%d) location(%s:%d) copyNumber(%.3f)", i, chromosome, PostitionStart, Values.get(i)));
                }
            }
            else
            {
                for(int i = 0; i < Values.size(); ++i)
                {
                    PPL_LOGGER.debug(format("index(%d) location(%s:%d) copyNumber(%.3f)", i, chromosome, PostitionStart, Values.get(i)));
                }
            }
        }
    }

    private class SampleCopyNumberTask implements Callable<Void>
    {
        private final List<String> mSampleIds;

        public SampleCopyNumberTask()
        {
            mSampleIds = Lists.newArrayList();
        }

        public List<String> getSampleIds() { return mSampleIds; }

        @Override
        public Void call()
        {
            for(int i = 0; i < mSampleIds.size(); ++i)
            {
                String sampleId = mSampleIds.get(i);

                processSample(sampleId);

                if(i > 0 && (i % 10) == 0)
                {
                    PPL_LOGGER.debug("processed {} samples", i);
                }
            }

            PPL_LOGGER.info("tasks complete for {} samples", mSampleIds.size());

            return null;
        }

        private void processSample(final String sampleId)
        {
            String samplePurpleDir = mPurpleDataDir.contains("*") ? mPurpleDataDir.replaceAll("\\*", sampleId) : mPurpleDataDir;

            try
            {
                PurplePurity purityData = PurplePurity.read(PurplePurity.generateFilename(samplePurpleDir, sampleId));
                double ploidy = purityData.Ploidy;

                List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(
                        PurpleCopyNumberFile.generateFilenameForReading(samplePurpleDir, sampleId));

                Map<String,List<Double>> chrCopyNumberValues = Maps.newHashMap();

                List<Double> currentCnValues = null;
                String currentChromosome = "";
                int lastPosition = 1;

                for(PurpleCopyNumber copyNumber : copyNumbers)
                {
                    if(!copyNumber.chromosome().equals(currentChromosome))
                    {
                        currentChromosome = copyNumber.chromosome();
                        currentCnValues = Lists.newArrayList();
                        chrCopyNumberValues.put(currentChromosome, currentCnValues);
                        lastPosition = 1;
                    }

                    double adjustedCopyNumber = copyNumber.averageTumorCopyNumber() / ploidy;

                    while(lastPosition + WINDOW_SIZE < copyNumber.end())
                    {
                        currentCnValues.add(adjustedCopyNumber);
                        lastPosition += WINDOW_SIZE;
                    }
                }

                addSampleCopyNumberValues(sampleId, chrCopyNumberValues);
            }
            catch(IOException e)
            {
                PPL_LOGGER.error("sample({}) failed to load purple files form {}: {}", sampleId, samplePurpleDir, e.toString());
            }
        }
    }

    public static final String FLD_PERCENTILE = "Percentile";

    private void writePercentilesFile()
    {
        if(mChrWindowsMap.isEmpty())
        {
            PPL_LOGGER.error("empty chromosome copy-number map");
            System.exit(1);
        }

        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFilename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_CHROMOSOME);
            sj.add(FLD_POSITION_START);
            sj.add(FLD_PERCENTILE);

            writer.write(sj.toString());
            writer.newLine();

            String firstChromosome = mChrWindowsMap.keySet().iterator().next();
            RefGenomeVersion refGenomeVersion = firstChromosome.startsWith(CHR_PREFIX) ? V38 : V37;

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                if(chromosome.isAllosome())
                    continue;

                String chrStr = refGenomeVersion.versionedChromosome(chromosome.toString());
                List<WindowValues> windows = mChrWindowsMap.get(chrStr);

                double lastPercentile = -1;

                for(WindowValues windowValues : windows)
                {
                    double percentile = windowValues.percentile(TARGET_PERCENTILE_COPY_NUMBER);

                    if(mSampleIds.size() >= 50 && (percentile == 0 || percentile == 1))
                        windowValues.logValues(chrStr, 50);

                    if(mCollapseOutput)
                    {
                        if(lastPercentile >= 0 && abs(lastPercentile - percentile) < 0.005)
                            continue;
                    }

                    sj = new StringJoiner(TSV_DELIM);
                    sj.add(chrStr);
                    sj.add(String.valueOf(windowValues.PostitionStart));

                    sj.add(format("%.2f", percentile));

                    writer.write(sj.toString());
                    writer.newLine();

                    lastPercentile = percentile;
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to write cohort copy number percentiles file({}): {}", mOutputFilename, e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        addSampleIdFile(configBuilder, true);

        configBuilder.addPath(PURPLE_DIR_CFG, true, "Directory pattern for sample purple directory");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Cohort copy number window percentiles");
        configBuilder.addFlag(COLLAPSE_OUTPUT, "Collapse positions with matching percentiles");

        addLoggingOptions(configBuilder);
        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        CopyNumberWindowPercentiles cnWindowPercentiles = new CopyNumberWindowPercentiles(configBuilder);
        cnWindowPercentiles.run();
    }
}