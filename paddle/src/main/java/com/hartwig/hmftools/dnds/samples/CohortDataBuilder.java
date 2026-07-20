package com.hartwig.hmftools.dnds.samples;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.dnds.DndsCommon.APP_NAME;
import static com.hartwig.hmftools.dnds.DndsCommon.DN_LOGGER;
import static com.hartwig.hmftools.dnds.SampleMutationalLoad.cohortSampleMutationalLoadFilename;
import static com.hartwig.hmftools.dnds.SampleMutationalLoad.initialiseWriter;

import java.io.BufferedWriter;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.dnds.SomaticVariant;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class CohortDataBuilder
{
    private final List<String> mSampleIds;
    private final String mPurpleDir;
    private final String mOutputDir;
    private final int mThreads;
    private final DatabaseAccess mDbAccess;
    private final SampleDataLoader mSampleDataLoader;

    private BufferedWriter mMutLoadWriter;
    private BufferedWriter mVariantsWriter;
    private final AtomicInteger mProcessedCount = new AtomicInteger();

    public CohortDataBuilder(final ConfigBuilder configBuilder)
    {
        mSampleIds = loadSampleIdsFile(configBuilder);
        mPurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        mOutputDir = parseOutputDir(configBuilder);
        mThreads = parseThreads(configBuilder);

        mDbAccess = DatabaseAccess.createDatabaseAccess(configBuilder);

        mSampleDataLoader = new SampleDataLoader(mPurpleDir, mDbAccess, configBuilder.getValue(TARGET_REGIONS_BED));

        mMutLoadWriter = null;
        mVariantsWriter = null;

    }

    public void run()
    {
        mVariantsWriter = SomaticVariant.initialiseWriter(SomaticVariant.cohortDndsVariantsFilename(mOutputDir));
        mMutLoadWriter = initialiseWriter(cohortSampleMutationalLoadFilename(mOutputDir));

        DN_LOGGER.info("retrieving variants and mutation load for {} samples", mSampleIds.size());

        List<SampleTask> sampleTasks = Lists.newArrayList();

        if(mThreads > 1)
        {
            for(int taskId = 0; taskId < min(mSampleIds.size(), mThreads); ++taskId)
            {
                SampleTask sampleTask = new SampleTask(taskId, mSampleIds.size(), mThreads, mProcessedCount,
                        mSampleDataLoader, mMutLoadWriter, mVariantsWriter);

                sampleTasks.add(sampleTask);
            }

            int taskIndex = 0;
            for(String sampleId : mSampleIds)
            {
                if(taskIndex >= sampleTasks.size())
                    taskIndex = 0;

                sampleTasks.get(taskIndex).addSample(sampleId);

                ++taskIndex;
            }

            final List<Callable<Void>> callableList = sampleTasks.stream().collect(Collectors.toList());

            if(!TaskExecutor.executeTasks(callableList, mThreads))
                System.exit(1);
        }
        else
        {
            SampleTask sampleTask = new SampleTask(0, mSampleIds.size(), mThreads, mProcessedCount,
                    mSampleDataLoader, mMutLoadWriter, mVariantsWriter);

            sampleTask.addSamples(mSampleIds);
            sampleTasks.add(sampleTask);
            sampleTask.call();
        }

        closeBufferedWriter(mMutLoadWriter);
        closeBufferedWriter(mVariantsWriter);

        DN_LOGGER.info("DNDS sample data building complete");
    }

    public static void main(final String... args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        CohortDataBuilder cohortDataBuilder = new CohortDataBuilder(configBuilder);
        cohortDataBuilder.run();
    }

    private static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(SAMPLE_ID_FILE, true, SAMPLE_ID_FILE_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(TARGET_REGIONS_BED, true, TARGET_REGIONS_BED_DESC);
        DatabaseAccess.addDatabaseCmdLineArgs(configBuilder, false);
        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }
}
