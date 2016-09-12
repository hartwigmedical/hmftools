package com.hartwig.hmftools.healthcheckeranalyser;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class HealthCheckerAnalysisApplication {

    private static final String REPORT_OUTPUT_PATH_ARGS_DESC = "The path where reports are written to.";
    private static final String REPORT_OUTPUT_PATH = "reportout";

    public static void main(final String... args) throws ParseException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String reportOutputPath = cmd.getOptionValue(REPORT_OUTPUT_PATH);

        if (reportOutputPath == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Health-Checks-Analyser", options);
            System.exit(1);
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(REPORT_OUTPUT_PATH, true, REPORT_OUTPUT_PATH_ARGS_DESC);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
