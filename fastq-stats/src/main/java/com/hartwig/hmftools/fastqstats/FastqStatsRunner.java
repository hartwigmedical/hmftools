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
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Fastq-Stats", options);
        } else if (filePath != null) {
            final FastqTracker tracker = FastqStats.processFile(filePath);
            writeOutputToCSV("", tracker, csvOutPath);
        } else {
            final String flowcellName = getFlowcellName(dirPath);
            final File baseCallsDir = getBaseCallsDir(dirPath);
            final long startTime = System.currentTimeMillis();
            final FastqTracker tracker = FastqStats.processDir(baseCallsDir);
            LOGGER.info("Total time: " + (System.currentTimeMillis() - startTime) + "ms.");
            writeOutputToCSV(flowcellName, tracker, csvOutPath);
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

    private static void writeOutputToCSV(@NotNull String flowcellName, @NotNull FastqTracker tracker,
            @NotNull String csvOutPath) throws IOException {
        final BufferedWriter writer = new BufferedWriter(new FileWriter(csvOutPath, false));
        writer.write("Flowcell " + flowcellName + ", " + tracker.flowcell().yield() + ", "
                + tracker.flowcell().q30() * 100.0 / tracker.flowcell().yield() + "\n");
        for (final String laneName : tracker.lanes().keySet()) {
            final FastqData lane = tracker.lane(laneName);
            writer.write("Lane " + laneName + ", " + lane.yield() + ", " + lane.q30() * 100.0 / lane.yield() + "\n");
        }
        for (final String sampleName : tracker.samples().keySet()) {
            final FastqData sample = tracker.sample(sampleName);
            writer.write("Sample " + sampleName + ", " + sample.yield() + ", " + sample.q30() * 100.0 / sample.yield()
                    + "\n");
        }
        final FastqData undetermined = tracker.undetermined();
        writer.write("Sample Undetermined, " + undetermined.yield() + ", " + undetermined.q30() + "\n");
        writer.close();
        LOGGER.info("Written fastq qc data to " + csvOutPath);
    }

    /**
     * Looks at the directory name and tries to extract the flowcell name.
     * Assumes the directory name is composed of multiple parts separated by the '_' character and the flowcell name is
     * on the 4th position.
     * Example of a valid directory: 170101_TEST_000_FLOWCELLNAME.
     *
     * @param dirPath path to the directory
     * @return The flowcell name, or "Unknown" if the directory name has less than 4 parts
     */

    @NotNull
    static String getFlowcellName(@NotNull String dirPath) {
        final File dir = new File(dirPath);
        final String[] dirNameArray = dir.getName().split("_");
        if (dirNameArray.length >= 4) {
            return dirNameArray[3];
        } else {
            LOGGER.warn("Could not get flowcell name from " + dir.getName());
            return "Unknown";
        }
    }

    @NotNull
    static File getBaseCallsDir(@NotNull String dirPath) throws IOException {
        final File baseCallsDir = new File(
                dirPath + File.separator + "Data" + File.separator + "Intensities" + File.separator + "BaseCalls");
        if (!baseCallsDir.exists()) {
            throw new IOException("dir " + baseCallsDir + " does not exist.");
        } else if (!baseCallsDir.isDirectory()) {
            throw new IOException(baseCallsDir + " is not a directory.");
        }
        return baseCallsDir;
    }

}
