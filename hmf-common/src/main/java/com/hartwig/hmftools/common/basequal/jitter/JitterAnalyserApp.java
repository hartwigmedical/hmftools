package com.hartwig.hmftools.common.basequal.jitter;

import static java.time.format.DateTimeFormatter.ISO_ZONED_DATE_TIME;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;

import java.io.IOException;
import java.time.Duration;
import java.time.Instant;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class JitterAnalyserApp
{
    public static final Logger sLogger = LogManager.getLogger(JitterAnalyserApp.class);

    private final JitterAnalyserConfig mConfig;

    public JitterAnalyserApp(final ConfigBuilder configBuilder) throws ParseException
    {
        mConfig = new JitterAnalyserConfig(configBuilder);
    }

    public int run(final String... args) throws InterruptedException, IOException
    {
        Instant start = Instant.now();

        VersionInfo versionInfo = new VersionInfo("errorprofile.version");

        sLogger.info("ErrorProfile version: {}", versionInfo.version());

        sLogger.debug("build timestamp: {}, run args: {}",
                versionInfo.buildTime().format(ISO_ZONED_DATE_TIME), String.join(" ", args));

        if(!mConfig.isValid())
        {
            sLogger.error(" invalid config, exiting");
            return 1;
        }

        JitterAnalyser jitterAnalyser = new JitterAnalyser(mConfig, sLogger);
        jitterAnalyser.processBam();
        jitterAnalyser.writeAnalysisOutput();

        Instant finish = Instant.now();
        long seconds = Duration.between(start, finish).getSeconds();
        sLogger.info("run complete. Time taken: {}m {}s", seconds / 60, seconds % 60);

        return 0;
    }

    public static void main(final String... args) throws InterruptedException, IOException, ParseException
    {
        ConfigBuilder configBuilder = new ConfigBuilder("ErrorProfile");
        JitterAnalyserConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        // set all thread exception handler
        // if we do not do this, exception thrown in other threads will not be handled and results
        // in the program hanging
        Thread.setDefaultUncaughtExceptionHandler((Thread t, Throwable e) ->
        {
            sLogger.error("[{}]: uncaught exception: {}", t, e);
            e.printStackTrace(System.err);
            System.exit(1);
        });

        System.exit(new JitterAnalyserApp(configBuilder).run(args));
    }
}
