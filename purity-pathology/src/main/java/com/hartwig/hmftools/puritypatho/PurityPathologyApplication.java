package com.hartwig.hmftools.puritypatho;

import java.io.IOException;

import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.puritypatho.variants.ReadingData;

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
    private static final String COUNT_SET = "count_set";
    private static final Logger LOGGER = LogManager.getLogger(PurityPathologyApplication.class);

    public static void main(final String... args) throws ParseException, IOException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String runsFolderPath = cmd.getOptionValue(RUNS_DIR);
        final String countSet = cmd.getOptionValue(COUNT_SET);
        LOGGER.info("runsFolderPath: " + runsFolderPath);
        final RunContext runContext = ProductionRunContextFactory.fromRunDirectory(runsFolderPath);
        final String tumorSample = runContext.tumorSample();
        LOGGER.info("tumorSample: " + tumorSample);
        LOGGER.info("countSet: " + countSet);
        ReadingData.readingFiles(runsFolderPath, tumorSample, countSet);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(RUNS_DIR, true, "Path towards the folder containing patient runs.");
        options.addOption(COUNT_SET, true, "Count of set");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}

