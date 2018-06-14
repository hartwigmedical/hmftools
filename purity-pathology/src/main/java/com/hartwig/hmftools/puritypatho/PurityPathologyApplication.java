package com.hartwig.hmftools.puritypatho;

import java.io.IOException;
import java.sql.SQLException;
import java.util.concurrent.ExecutionException;

import com.hartwig.hmftools.common.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PurityPathologyApplication {
    private static final String RUNS_DIR = "runs_dir";

    private static final Logger LOGGER = LogManager.getLogger(PurityPathologyApplication.class);

    public static void main(final String... args) throws ParseException, IOException, SQLException, ExecutionException, InterruptedException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String runsFolderPath = cmd.getOptionValue(RUNS_DIR);
        LOGGER.info("Set name: " + runsFolderPath);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(RUNS_DIR, true, "Path towards the folder containing patient runs.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}

