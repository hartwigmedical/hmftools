package com.hartwig.hmftools.purple.tools;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.purple.GermlineStatus.DIPLOID;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.segment.SegmentFile;

import org.jetbrains.annotations.NotNull;

public class SegmentAnalyser
{
    private final List<String> mSampleIds;
    private final String mPurpleDir;
    private final int mThreads;

    private final BufferedWriter mWriter;

    public SegmentAnalyser(final ConfigBuilder configBuilder)
    {
        mSampleIds = loadSampleIdsFile(configBuilder);
        mPurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        mThreads = parseThreads(configBuilder);

        mWriter = initialiseWriter(parseOutputDir(configBuilder));
    }

    public void run()
    {
        if(mSampleIds.isEmpty())
        {
            PPL_LOGGER.error("missing sampleIds, exiting");
            System.exit(1);
        }

        PPL_LOGGER.info("running Purple segment analyser for {} samples", mSampleIds.size());

        List<SampleTask> sampleTasks = Lists.newArrayList();

        for(int i = 0; i < min(mSampleIds.size(), mThreads); ++i)
        {
            sampleTasks.add(new SampleTask(i));
        }

        int taskIndex = 0;
        for(String sampleId : mSampleIds)
        {
            if(taskIndex >= sampleTasks.size())
                taskIndex = 0;

            sampleTasks.get(taskIndex).getSampleIds().add(sampleId);

            ++taskIndex;
        }

        final List<Callable> callableList = sampleTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mThreads);

        closeBufferedWriter(mWriter);

        PPL_LOGGER.info("Purple segment analysis complete");
    }

    private class SampleTask implements Callable
    {
        private final int mTaskId;
        private final List<String> mSampleIds;

        public SampleTask(int taskId)
        {
            mTaskId = taskId;
            mSampleIds = Lists.newArrayList();
        }

        public List<String> getSampleIds() { return mSampleIds; }

        @Override
        public Long call()
        {
            for(int i = 0; i < mSampleIds.size(); ++i)
            {
                String sampleId = mSampleIds.get(i);

                processSample(sampleId);

                if(i > 0 && (i % 100) == 0)
                {
                    PPL_LOGGER.debug("{}: processed {} samples", mTaskId, i);
                }
            }

            PPL_LOGGER.info("{}: tasks complete for {} samples", mTaskId, mSampleIds.size());

            return (long)0;
        }

        private void processSample(final String sampleId)
        {
            String samplePurpleDir = convertWildcardSamplePath(mPurpleDir, sampleId);

            try
            {
                List<ObservedRegion> fittedRegions = SegmentFile.read(SegmentFile.generateFilename(samplePurpleDir, sampleId));

                double minRatio = 0;
                double maxRatio = 0;

                for(ObservedRegion region : fittedRegions)
                {
                    if(ignoreRegion(region, true))
                        continue;

                    double normalRatio = region.observedNormalRatio() / region.unnormalisedObservedNormalRatio();

                    maxRatio = max(maxRatio, normalRatio);
                    minRatio = minRatio > 0 ? min(minRatio, normalRatio) : normalRatio;
                }

                double logMaxThreshold = maxRatio >= 1.3 ? 0.9 * maxRatio : maxRatio;
                double logMinThreshold = minRatio <= 0.8 ? 1.1 * minRatio : 0;

                for(ObservedRegion region : fittedRegions)
                {
                    if(ignoreRegion(region, false))
                        continue;

                    double normalRatio = region.observedNormalRatio() / region.unnormalisedObservedNormalRatio();

                    if(normalRatio >= logMaxThreshold || normalRatio <= logMinThreshold)
                        writeSegmentData(sampleId, region);
                }
            }
            catch(IOException e)
            {
                PPL_LOGGER.error("sample({}) failed to load Purple segment file form {}: {}", sampleId, samplePurpleDir, e.toString());
            }
        }
    }

    private static boolean ignoreRegion(final ObservedRegion region, boolean onlyDiploid)
    {
        if(HumanChromosome.fromString(region.chromosome()).isAllosome())
            return true;

        if(onlyDiploid)
        {
            if(region.germlineStatus() != DIPLOID)
                return true;
        }
        else
        {
            if(region.germlineStatus() != DIPLOID && region.germlineStatus() != HET_DELETION && region.germlineStatus() != HOM_DELETION)
                return true;
        }

        return region.unnormalisedObservedNormalRatio() == 0;
    }

    private BufferedWriter initialiseWriter(final String outputDir)
    {
        try
        {
            String fileName = outputDir + "purple_segment_analysis.csv";

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("SampleId\tChromosome\tRegionStart\tRegionEnd\tGermlineStatus\tBafCount\tDepthWindows");
            writer.write("\tObsTumorRatio\tObsNormalRatio\tObsUnnormalisedNormalRatio\tRegionRefNormalisedCN");


            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to initialise output file output: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeSegmentData(final String sampleId, final ObservedRegion region)
    {
        try
        {
            mWriter.write(String.format("%s\t%s\t%d\t%d\t%s",
                    sampleId, region.chromosome(), region.start(), region.end(), region.germlineStatus()));

            mWriter.write(String.format("\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f",
                    region.bafCount(), region.depthWindowCount(), region.observedTumorRatio(),
                    region.observedNormalRatio(), region.unnormalisedObservedNormalRatio(), region.refNormalisedCopyNumber()));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to write germline gene deletion file output: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        addSampleIdFile(configBuilder, true);

        configBuilder.addPath(PURPLE_DIR_CFG, true, "Directory pattern for sample purple directory");

        addLoggingOptions(configBuilder);
        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        SegmentAnalyser segmentAnalyser = new SegmentAnalyser(configBuilder);
        segmentAnalyser.run();
    }
}