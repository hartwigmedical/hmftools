package com.hartwig.hmftools.fastqstats;

import java.io.File;
import java.io.IOException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public final class FastqStatsRunner {
    private static final Logger LOGGER = LogManager.getLogger(FastqStatsRunner.class);
    private static final String FASTQ_FILE = "file";
    private static final String FASTQ_ROOT_DIR = "dir";

    public static void main(String[] args) throws ParseException, IOException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String fileName = cmd.getOptionValue(FASTQ_FILE);
        final String dirName = cmd.getOptionValue(FASTQ_ROOT_DIR);

        if (fileName == null && dirName == null) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Fastq-Stats", options);
        }
        else if(fileName != null){
            FastqTracker tr = FastqStats.processFile(fileName);
            System.out.println(tr);
        }
        else {
            File dir = new File(dirName);
            if(dir.isDirectory()) {
                FastqTracker tr = FastqStats.processDir(dir);
                System.out.println(tr);
            } else {
                if(!dir.exists()){
                    LOGGER.warn("dir " + dir + " does not exist.");
                }
                HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp("Fastq-Stats", options);
            }
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(FASTQ_FILE, true, "Path towards the original fastq file.");
        options.addOption(FASTQ_ROOT_DIR, true, "Path towards the root fastq dir.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
