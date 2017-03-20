package com.hartwig.hmftools.fastqstats;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
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

public final class FastqStatsRunner {
    private static final Logger LOGGER = LogManager.getLogger(FastqStatsRunner.class);
    private static final String FASTQ_FILE = "file";
    private static final String FASTQ_ROOT_DIR = "dir";
    private static final String CSV_OUT_DIR = "out";

    public static void main(String[] args) throws ParseException, IOException, InterruptedException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String filePath = cmd.getOptionValue(FASTQ_FILE);
        final String dirPath = cmd.getOptionValue(FASTQ_ROOT_DIR);
        final String csvOutPath = cmd.getOptionValue(CSV_OUT_DIR);

        if ((filePath == null && dirPath == null) || csvOutPath == null) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Fastq-Stats", options);
        } else if (filePath != null) {
            FastqTracker tr = FastqStats.processFile(filePath);
            writeOutputToCSV("", tr, csvOutPath);
        } else {
            File dir = new File(dirPath);
            File baseCallsDir = new File(dirPath + "/Data/Intensities/BaseCalls/");
            String[] dirNameArray = dir.getName().split("_");
            String flowcellName;
            if (dirNameArray.length >= 4) {
                flowcellName = dirNameArray[3];
            } else {
                LOGGER.warn("Could not get flowcell name from " + dir.getName());
                flowcellName = "Unknown";
            }
            if (baseCallsDir.isDirectory()) {
                final long startTime = System.currentTimeMillis();
                FastqTracker tr = FastqStats.processDir(baseCallsDir);
                LOGGER.info("Total time: " + (System.currentTimeMillis() - startTime) + "ms.");
                writeOutputToCSV(flowcellName, tr, csvOutPath);
            } else {
                if (!baseCallsDir.exists()) {
                    LOGGER.warn("dir " + baseCallsDir + " does not exist.");
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
        options.addOption(FASTQ_ROOT_DIR, true, "Path towards the flowcell dir.");
        options.addOption(CSV_OUT_DIR, true, "Path towards the csv output file.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    public static void writeOutputToCSV(@NotNull String flowcellName, @NotNull FastqTracker tracker,
            @NotNull String csvOutPath) throws IOException {
        final BufferedWriter writer = new BufferedWriter(new FileWriter(csvOutPath, false));
        writer.write("Flowcell " + flowcellName + ", " + tracker.getFlowcellData().getYield() + ", "
                + tracker.getFlowcellData().getQ30() * 100.0 / tracker.getFlowcellData().getYield() + "\n");
        for (String laneName : tracker.getLanes().keySet()) {
            FastqData lane = tracker.getLaneData(laneName);
            writer.write("Lane " + laneName + ", " + lane.getYield() + ", " + lane.getQ30() * 100.0 / lane.getYield()
                    + "\n");
        }
        for (String sampleName : tracker.getSamples().keySet()) {
            FastqData sample = tracker.getSampleData(sampleName);
            writer.write("Sample " + sampleName + ", " + sample.getYield() + ", "
                    + sample.getQ30() * 100.0 / sample.getYield() + "\n");
        }
        FastqData undetermined = tracker.getUndeterminedData();
        writer.write("Sample Undetermined, " + undetermined.getYield() + ", " + undetermined.getQ30() + "\n");
        writer.close();
        LOGGER.info("Written fastq qc data to " + csvOutPath);
    }
}
