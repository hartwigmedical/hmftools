package com.hartwig.hmftools.puritypatho;

import java.io.IOException;

import com.hartwig.hmftools.puritypatho.coverages.ReadingFastqFiles;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CoveragesFastqApplication {
    private static final String FASTQ_DIR = "fastq_dir";
    private static final Logger LOGGER = LogManager.getLogger(CoveragesFastqApplication.class);

    public static void main(final String... args) throws ParseException, IOException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String fastqDir = cmd.getOptionValue(FASTQ_DIR);
        LOGGER.info("fastqDir application: " + fastqDir);
        ReadingFastqFiles.readingFiles(fastqDir);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(FASTQ_DIR, true, "Path towards the folder containing fastq files.");
        return options;
    }
    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
