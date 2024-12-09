package com.hartwig.hmftools.redux.jitter;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyser;
import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConfig;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.ParseException;

public class GenerateMsiJitterParams
{
    private final JitterAnalyserConfig mConfig;

    private final String mBamPath;
    private final int mThreads;
    private final int mPartitionSize;

    private static final String BAM_FILE = "bam";
    private static final String PARTITION_SIZE = "partition_size";

    public static final int DEFAULT_PARTITION_SIZE = 1_000_000;

    public GenerateMsiJitterParams(final ConfigBuilder configBuilder) throws ParseException
    {
        mConfig = new JitterAnalyserConfig(configBuilder);

        mBamPath = configBuilder.getValue("bam");
        mPartitionSize = configBuilder.getInteger(PARTITION_SIZE);
        mThreads = parseThreads(configBuilder);
    }

    public void run() throws InterruptedException, IOException
    {
        long startTimeMs = System.currentTimeMillis();

        VersionInfo versionInfo = new VersionInfo("errorprofile.version");

        RD_LOGGER.info("ErrorProfile version: {}", versionInfo.version());

        if(!mConfig.isValid())
        {
            RD_LOGGER.error(" invalid config, exiting");
            System.exit(1);
        }

        JitterAnalyser jitterAnalyser = new JitterAnalyser(mConfig, RD_LOGGER, null);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("worker-%d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mThreads, namedThreadFactory);

        SampleBamProcessor sampleBamProcessor = new SampleBamProcessor(
                mConfig, mPartitionSize, jitterAnalyser.bamSlicerFilter(), jitterAnalyser.sampleReadProcessor());

        sampleBamProcessor.queryBam(mConfig, executorService, mBamPath);

        jitterAnalyser.writeAnalysisOutput();

        RD_LOGGER.info("Redux MSi jitter site analysis complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public static void main(final String... args) throws InterruptedException, IOException, ParseException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        JitterAnalyserConfig.registerConfig(configBuilder);

        configBuilder.addPath(BAM_FILE, true, "Path to bam file");
        configBuilder.addInteger(PARTITION_SIZE, "Size of the partitions processed by worker threads", DEFAULT_PARTITION_SIZE);

        addThreadOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);


        new GenerateMsiJitterParams(configBuilder).run();
    }
}
