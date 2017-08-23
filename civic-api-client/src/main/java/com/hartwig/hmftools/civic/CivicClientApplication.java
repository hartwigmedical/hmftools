package com.hartwig.hmftools.civic;

import com.hartwig.hmftools.civic.api.CivicApiWrapper;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CivicClientApplication {

    private static final Logger LOGGER = LogManager.getLogger(CivicClientApplication.class);

    public static void main(final String... args) throws Exception {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        CivicApiWrapper.getAllEvidenceItems().blockingSubscribe(LOGGER::info, LOGGER::error, () -> LOGGER.info("completed"));
        CivicApiWrapper.releaseResources();
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
