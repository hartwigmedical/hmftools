package com.hartwig.hmftools.dnds.builder;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.dnds.DndsCommon.APP_NAME;
import static com.hartwig.hmftools.dnds.DndsCommon.DN_LOGGER;
import static com.hartwig.hmftools.dnds.DndsCommon.SOMATIC_CACHE_DIR;
import static com.hartwig.hmftools.dnds.SampleMutationalLoad.cohortSampleMutationalLoadFilename;
import static com.hartwig.hmftools.dnds.SampleMutationalLoad.initialiseWriter;
import static com.hartwig.hmftools.dnds.SampleMutationalLoad.loadCohortSampleMutationalLoads;

import java.io.BufferedWriter;
import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.dnds.SampleMutationalLoad;
import com.hartwig.hmftools.dnds.SomaticVariant;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class DndsDataBuilder
{
    private final List<String> mSampleIds;
    private final String mPurpleDir;
    private final String mOutputDir;
    private final int mThreads;
    private final DatabaseAccess mDbAccess;
    private final SampleDataLoader mSampleDataLoader;

    private BufferedWriter mMutLoadWriter;
    private final String mSomaticsDir;

    public DndsDataBuilder(final ConfigBuilder configBuilder)
    {
        mSampleIds = loadSampleIdsFile(configBuilder);
        mPurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        mThreads = parseThreads(configBuilder);
        mOutputDir = parseOutputDir(configBuilder);
        mSomaticsDir = mOutputDir + SOMATIC_CACHE_DIR + File.separator;
        mMutLoadWriter = null;

        mDbAccess = DatabaseAccess.createDatabaseAccess(configBuilder);
        mSampleDataLoader = new SampleDataLoader(mPurpleDir, mDbAccess);
    }

    public void run()
    {
        DN_LOGGER.info("DNDS sample data building for {} samples", mSampleIds.size());

        // load exonic somatics either from file or the DB
        String somaticsDir = mOutputDir + SOMATIC_CACHE_DIR + File.separator;
        new File(somaticsDir).mkdirs();

        // check for existing sample data
        String cohortSampleMutLoadFilename = cohortSampleMutationalLoadFilename(mOutputDir);

        Map<String, SampleMutationalLoad> sampleMutationalLoadMap = loadCohortSampleMutationalLoads(cohortSampleMutLoadFilename, false);

        List<String> missingSampleIds = mSampleIds.stream().filter(x -> !sampleMutationalLoadMap.containsKey(x)).collect(Collectors.toList());

        if(missingSampleIds.isEmpty())
        {
            DN_LOGGER.info("all {} samples have mutational load data", mSampleIds.size());
            return;
        }

        mMutLoadWriter = initialiseWriter(cohortSampleMutLoadFilename);

        DN_LOGGER.info("retrieving mutation load and variant data for {} samples", missingSampleIds.size());

        List<SampleTask> sampleTasks = Lists.newArrayList();

        if(mThreads > 1)
        {
            for(int i = 0; i < min(missingSampleIds.size(), mThreads); ++i)
            {
                sampleTasks.add(new SampleTask(i));
            }

            int taskIndex = 0;
            for(String sampleId : missingSampleIds)
            {
                if(taskIndex >= sampleTasks.size())
                    taskIndex = 0;

                sampleTasks.get(taskIndex).addSample(sampleId);

                ++taskIndex;
            }

            final List<Callable> callableList = sampleTasks.stream().collect(Collectors.toList());

            if(!TaskExecutor.executeTasks(callableList, mThreads))
                System.exit(1);
        }
        else
        {
            SampleTask sampleTask = new SampleTask(0);
            sampleTask.addSamples(missingSampleIds);
            sampleTasks.add(sampleTask);
            sampleTask.call();
        }

        closeBufferedWriter(mMutLoadWriter);

        DN_LOGGER.info("DNDS sample data building complete");
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

        public void addSample(final String sampleId) { mSampleIds.add(sampleId); }
        public void addSamples(final List<String> sampleIds) { mSampleIds.addAll(sampleIds); }

        @Override
        public Long call()
        {
            for(int i = 0; i < mSampleIds.size(); ++i)
            {
                String sampleId = mSampleIds.get(i);

                processSample(sampleId);
                if(i > 0 && (i % 10) == 0)
                {
                    DN_LOGGER.info("{}: processed {} samples", mTaskId, i);
                }
            }

            if(mThreads > 1)
            {
                DN_LOGGER.info("{}: tasks complete for {} samples", mTaskId, mSampleIds.size());
            }

            return (long)0;
        }

        private void processSample(final String sampleId)
        {
            List<SomaticVariant> variants = mSampleDataLoader.loadVariants(sampleId);
            DN_LOGGER.debug("sample({}) loaded {} variants", sampleId, variants.size());

            SampleMutationalLoad sampleMutationalLoad = mSampleDataLoader.calcSampleMutationalLoad(sampleId);
            SampleMutationalLoad.writeSampleMutationalLoad(mMutLoadWriter, sampleId, sampleMutationalLoad);
            SomaticVariant.writeVariants(mSomaticsDir, sampleId, variants);
        }
    }

    public static void main(final String... args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        DndsDataBuilder dndsDataBuilder = new DndsDataBuilder(configBuilder);
        dndsDataBuilder.run();
    }

    private static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(SAMPLE_ID_FILE, true, SAMPLE_ID_FILE_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        DatabaseAccess.addDatabaseCmdLineArgs(configBuilder, false);
        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }
}
