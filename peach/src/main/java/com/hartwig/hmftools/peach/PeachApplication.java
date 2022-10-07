package com.hartwig.hmftools.peach;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;

import org.apache.commons.cli.*;
import org.jetbrains.annotations.NotNull;

import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;

public class PeachApplication {
    private final PeachConfig mConfig;
    public PeachApplication(final PeachConfig config)
    {
        mConfig = config;
    }

    public void run(){
        PCH_LOGGER.info("running PEACH");
    }

    public static void main(String[] args) {
        final Options options = PeachConfig.createOptions();

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            PeachConfig config = new PeachConfig(cmd);
            PeachApplication peachApplication = new PeachApplication(config);
            peachApplication.run();
        }
        catch(ParseException e)
        {
            PCH_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("PeachApplication", options);
            System.exit(1);
        }
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}