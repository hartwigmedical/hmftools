package com.hartwig.hmftools.dnds.builder;

import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.dnds.DndsCommon.DN_LOGGER;
import static com.hartwig.hmftools.dnds.DndsCommon.SOMATIC_CACHE_DIR;
import static com.hartwig.hmftools.dnds.DndsCommon.logVersion;
import static com.hartwig.hmftools.dnds.SampleMutationalLoad.cohortSampleMutationalLoadFilename;
import static com.hartwig.hmftools.dnds.SampleMutationalLoad.initialiseWriter;
import static com.hartwig.hmftools.dnds.SampleMutationalLoad.loadCohortSampleMutationalLoads;

import java.io.BufferedWriter;
import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

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

    public DndsDataBuilder(final ConfigBuilder configBuilder)
    {
        mSampleIds = loadSampleIdsFile(configBuilder);
        mPurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        mThreads = parseThreads(configBuilder);
        mOutputDir = parseOutputDir(configBuilder);

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

        BufferedWriter mutLoadWriter = null;

        List<String> missingSampleIds = mSampleIds.stream().filter(x -> !sampleMutationalLoadMap.containsKey(x)).collect(Collectors.toList());

        if(!missingSampleIds.isEmpty())
        {
            mutLoadWriter = initialiseWriter(cohortSampleMutLoadFilename);

            DN_LOGGER.info("retrieving mutation load and variant data for {} samples", missingSampleIds.size());

            for(String sampleId : missingSampleIds)
            {
                processSample(sampleId, somaticsDir, mutLoadWriter);
            }

            closeBufferedWriter(mutLoadWriter);
        }

        DN_LOGGER.info("DNDS sample data building complete");
    }

    private void processSample(final String sampleId, final String somaticsDir, final BufferedWriter mutLoadWriter)
    {
        List<SomaticVariant> variants = mSampleDataLoader.loadVariants(sampleId);
        DN_LOGGER.debug("sample({}) loaded {} variants", sampleId, variants.size());

        SampleMutationalLoad sampleMutationalLoad = mSampleDataLoader.calcSampleMutationalLoad(sampleId);
        SampleMutationalLoad.writeSampleMutationalLoad(mutLoadWriter, sampleId, sampleMutationalLoad);
        SomaticVariant.writeVariants(somaticsDir, sampleId, variants);
    }

    public static void main(final String... args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);
        setLogLevel(configBuilder);

        logVersion();

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
