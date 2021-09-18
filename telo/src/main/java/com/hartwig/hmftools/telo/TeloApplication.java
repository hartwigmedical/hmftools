package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;

import java.time.Instant;
import java.time.Duration;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class TeloApplication
{
    private final TeloConfig mConfig;

    public TeloApplication(final Options options, final String... args) throws ParseException
    {
        VersionInfo versionInfo = new VersionInfo("telo.version");
        TE_LOGGER.info("Telo version: {}", versionInfo.version());

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        mConfig = new TeloConfig(cmd);
    }

    private void run()
    {
        if(!mConfig.isValid())
        {
            TE_LOGGER.error(" invalid config, exiting");
            System.exit(1);
        }

        TE_LOGGER.info("starting telomeric analysis");
        Instant start = Instant.now();

        try
        {
            BamProcessor.processBam(mConfig);
        }
        catch (InterruptedException e)
        {
            TE_LOGGER.error("Telo run interrupted, exiting");
            System.exit(1);
        }

        Instant finish = Instant.now();
        long seconds = Duration.between(start, finish).getSeconds();
        TE_LOGGER.info("Telo run complete, time taken: {}m {}s",  seconds / 60, seconds % 60);
    }

    public static void main(final String... args)
    {
        final Options options = TeloConfig.createOptions();
        TeloApplication application = null;

        try
        {
            application = new TeloApplication(options, args);
        }
        catch(ParseException e)
        {
            TE_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("TeloApplication", options);
            System.exit(1);
        }

        application.run();
    }

    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
