package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.ExecutionException;

import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class TeloApplication
{
    private final TeloConfig mConfig;

    public TeloApplication(final Options options, final String... args) throws ParseException, IOException
    {
        VersionInfo versionInfo = new VersionInfo("cobalt.version");
        TE_LOGGER.info("Telo version: {}", versionInfo.version());

        final CommandLine cmd = createCommandLine(args, options);
        mConfig = new TeloConfig(cmd);
    }

    private void run()
    {
        if(!mConfig.isValid())
        {
            TE_LOGGER.error(" invalid config, exiting");
            System.exit(1);
        }

        TE_LOGGER.info("starting Telo");

        BamReader bamReader = new BamReader(mConfig);

        bamReader.findTelomereContent();

        TE_LOGGER.info("telo run complete");
    }

    public static void main(final String... args) throws IOException
    {
        final Options options = TeloConfig.createOptions();

        try
        {
            TeloApplication application = new TeloApplication(options, args);
            application.run();
        }
        catch(ParseException e)
        {
            TE_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("CountBamLinesApplication", options);
            System.exit(1);
        }
    }

    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
