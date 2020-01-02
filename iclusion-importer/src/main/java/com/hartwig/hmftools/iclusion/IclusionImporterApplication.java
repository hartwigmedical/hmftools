package com.hartwig.hmftools.iclusion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class IclusionImporterApplication {
    private static final Logger LOGGER = LogManager.getLogger(IclusionImporterApplication.class);

    public static void main(@NotNull final String[] args) throws ParseException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();

        return options;
    }

}
