package com.hartwig.hmftools.dnds.calcs;

import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.dnds.DndsCommon.APP_NAME;
import static com.hartwig.hmftools.dnds.DndsCommon.DN_LOGGER;
import static com.hartwig.hmftools.dnds.DndsCommon.SOMATIC_CACHE_DIR;
import static com.hartwig.hmftools.dnds.SampleMutationalLoad.cohortSampleMutationalLoadFilename;
import static com.hartwig.hmftools.dnds.SampleMutationalLoad.loadCohortSampleMutationalLoads;
import static com.hartwig.hmftools.dnds.SomaticVariant.cohortDndsVariantsFilename;
import static com.hartwig.hmftools.dnds.calcs.CohortMutationalLoad.fromSampleMap;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.dnds.SampleMutationalLoad;

public class DndsFileBuilder
{
    private final List<String> mSampleIds;
    private final String mOutputDir;
    private final int mThreads;

    public DndsFileBuilder(final ConfigBuilder configBuilder)
    {
        mSampleIds = loadSampleIdsFile(configBuilder);
        mThreads = parseThreads(configBuilder);
        mOutputDir = parseOutputDir(configBuilder);
    }

    public void run()
    {
        DN_LOGGER.info("DNDS file builder from {} samples", mSampleIds.size());

        // load exonic somatics either from file or the DB
        String somaticsDir = mOutputDir + SOMATIC_CACHE_DIR + File.separator;

        if(!Files.exists(Paths.get(somaticsDir)))
        {
            DN_LOGGER.error("missing somatics directory", somaticsDir);
            System.exit(1);
        }

        Map<String,SampleMutationalLoad> sampleMutationalLoadMap = loadCohortSampleMutationalLoads(
                cohortSampleMutationalLoadFilename(mOutputDir), true);

        CohortMutationalLoad cohortMutationalLoad = fromSampleMap(sampleMutationalLoadMap);

        Map<String,DndsCvGene> geneDndsCvmap = DndsCvGene.load(DndsCvGene.geneDndsCvFilename(mOutputDir));

        List<DndsMutation> cohortMutations = DndsMutation.readVariants(cohortDndsVariantsFilename(mOutputDir));


        DN_LOGGER.info("DNDS file building complete");
    }

    public static void main(final String... args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        DndsFileBuilder dndsFileBuilder = new DndsFileBuilder(configBuilder);
        dndsFileBuilder.run();
    }

    private static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(SAMPLE_ID_FILE, true, SAMPLE_ID_FILE_DESC);
        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }
}
