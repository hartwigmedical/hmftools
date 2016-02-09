package com.hartwig.hmftools.sullivan;

import org.apache.commons.cli.*;
import org.jetbrains.annotations.NotNull;

public class SullivanRunner {

    private static final String ORIGINAL_FASTQ_PATH_VAR = "origpath";
    private static final String RECREATED_FASTQ_PATH_VAR = "recpath";
    private static final String DIRECTORY_MODE_VAR = "d";
    private static final String NUM_RECORDS_VAR = "n";

    private static final int DEFAULT_NUM_RECORDS = 100000;

    public static void main(String[] args) throws ParseException {
        Options options = createOptions();
        CommandLine cmd = createCommandLine(args, options);

        String originalPath = cmd.getOptionValue(ORIGINAL_FASTQ_PATH_VAR);
        String recreatedPath = cmd.getOptionValue(RECREATED_FASTQ_PATH_VAR);

        if (originalPath == null || recreatedPath == null) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("sullivan", options);
        } else {
            int numRecords =
                    Integer.parseInt(cmd.getOptionValue(NUM_RECORDS_VAR, Integer.toString(DEFAULT_NUM_RECORDS)));
            boolean isDirectoryMode = cmd.getOptionValue(DIRECTORY_MODE_VAR) != null;

            SullivanAlgo.runSullivanAlgo(originalPath, recreatedPath, isDirectoryMode, numRecords);
        }
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options)
            throws ParseException {
        CommandLineParser parser = new DefaultParser();
        return parser.parse(createOptions(), args);
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(ORIGINAL_FASTQ_PATH_VAR, true, "Path towards the original fastq file(s)");
        options.addOption(RECREATED_FASTQ_PATH_VAR, true, "Path towards the recreated fastq file(s)");
        options.addOption(NUM_RECORDS_VAR, "num-records", true, "Number of records from original file(s) to validate");
        options.addOption(DIRECTORY_MODE_VAR, "directory-mode", false, "If set, run sullivan on every file in the provided paths");
        return options;
    }
}
