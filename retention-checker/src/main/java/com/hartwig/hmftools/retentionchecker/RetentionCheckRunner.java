package com.hartwig.hmftools.retentionchecker;

import java.util.Collection;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class RetentionCheckRunner {

    private static final Logger LOGGER = LogManager.getLogger(RetentionCheckRunner.class);

    private static final String FASTQ_PATH1_VAR = "fastqpath1";
    private static final String FASTQ_PATH2_VAR = "fastqpath2";
    private static final String RECREATED_FASTQ_PATH_VAR = "recreated";
    private static final String NUM_RECORDS_VAR = "n";

    private static final int DEFAULT_NUM_RECORDS = 100000;

    public static void main(String[] args) throws ParseException {
        Options options = createOptions();
        CommandLine cmd = createCommandLine(args, options);

        String fastqPath1 = cmd.getOptionValue(FASTQ_PATH1_VAR);
        String recreatedPath = cmd.getOptionValue(RECREATED_FASTQ_PATH_VAR);

        if (fastqPath1 == null || recreatedPath == null) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("retention-checker", options);
        } else {
            int numRecords =
                    Integer.parseInt(cmd.getOptionValue(NUM_RECORDS_VAR, Integer.toString(DEFAULT_NUM_RECORDS)));

            Collection<String> fastqPaths = Lists.newArrayList(fastqPath1);
            String fastqPath2 = cmd.getOptionValue(FASTQ_PATH2_VAR);
            if (fastqPath2 != null) {
                fastqPaths.add(fastqPath2);
            }

            FileNameConverter converter = new DefaultNameConverter();
            RetentionCheckAlgo algo = new RetentionCheckAlgo(converter);

            boolean success = algo.runAlgo(fastqPaths, recreatedPath, numRecords);

            if (success) {
                LOGGER.info("RetentionCheck ran successfully!");
            } else {
                LOGGER.info("RetentionCheck not successful!");
            }
        }
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(FASTQ_PATH1_VAR, true, "Path towards the original fastq file(s)");
        options.addOption(FASTQ_PATH2_VAR, true, "Path towards the original fastq file(s)");
        options.addOption(RECREATED_FASTQ_PATH_VAR, true, "Path towards the recreated fastq file(s)");
        options.addOption(NUM_RECORDS_VAR, "num-records", true, "Number of records from original file(s) to validate");
        return options;
    }
}
