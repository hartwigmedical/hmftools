package com.hartwig.patientreporter;

import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PatientReporterApplication {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterApplication.class);

    private static final String RUN_DIRECTORY_ARGS_DESC = "A path towards a single rundir.";
    private static final String RUN_DIRECTORY = "rundir";

    @NotNull
    private final String runDirectory;

    public static void main(final String... args) throws ParseException, IOException {
        Options options = createOptions();
        CommandLine cmd = createCommandLine(options, args);

        String runDir = cmd.getOptionValue(RUN_DIRECTORY);

        if (runDir == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Health-Checks-Analyser", options);
            System.exit(1);
        }

        new PatientReporterApplication(runDir).printRunDirectory();
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(RUN_DIRECTORY, true, RUN_DIRECTORY_ARGS_DESC);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);

    }

    PatientReporterApplication(@NotNull final String runDirectory) {
        this.runDirectory = runDirectory;
    }

    void printRunDirectory() {
        LOGGER.info(runDirectory);
    }
}
