package com.hartwig.hmftools.sullivan;

import org.apache.commons.cli.*;
import org.jetbrains.annotations.NotNull;

import static com.hartwig.hmftools.sullivan.Logger.log;

public class SullivanRunner {

    private static final String ORIGINAL_FASTQ_PATH_VAR = "origpath";
    private static final String MERGE_ORIGINAL_FASTQ_PATH_VAR = "mergeorigpath";
    private static final String RECREATED_FASTQ_PATH_VAR = "recpath";
    private static final String DIRECTORY_MODE_VAR = "d";
    private static final String ANCIENT_FILE_NAME_MODE_VAR = "a";
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
            boolean isDirectoryMode = cmd.hasOption(DIRECTORY_MODE_VAR);
            String mergeOrigPath = cmd.getOptionValue(MERGE_ORIGINAL_FASTQ_PATH_VAR);
            boolean isAncientFileNameMode = cmd.hasOption(ANCIENT_FILE_NAME_MODE_VAR);

            FileNameConverter converter =
                    isAncientFileNameMode ? new AncientNameConverter() : new DefaultNameConverter();
            SullivanAlgo algo = new SullivanAlgo(converter);
            boolean success =
                    algo.runSullivanAlgo(originalPath, recreatedPath, mergeOrigPath, isDirectoryMode, numRecords);

            if (success) {
                log("Sullivan ran successfully!");
            } else {
                log("Sullivan not successful!");
            }
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
        options.addOption(MERGE_ORIGINAL_FASTQ_PATH_VAR, true, "Extra path towards the original fastq file(s), for merge runs");
        options.addOption(RECREATED_FASTQ_PATH_VAR, true, "Path towards the recreated fastq file(s)");
        options.addOption(NUM_RECORDS_VAR, "num-records", true, "Number of records from original file(s) to validate");
        options.addOption(DIRECTORY_MODE_VAR, "directory-mode", false, "If set, run sullivan on every file in the provided paths");
        options.addOption(ANCIENT_FILE_NAME_MODE_VAR, "ancient-file-name-mode", false, "If set, sullivan expects an ancient format for file names.");
        return options;
    }
}
