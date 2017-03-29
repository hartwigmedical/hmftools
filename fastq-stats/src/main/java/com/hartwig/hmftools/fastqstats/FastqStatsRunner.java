package com.hartwig.hmftools.fastqstats;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;

import com.google.common.annotations.VisibleForTesting;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.w3c.dom.Document;

public final class FastqStatsRunner {
    private static final Logger LOGGER = LogManager.getLogger(FastqStatsRunner.class);

    private static final String FASTQ_FILE = "file";
    private static final String FASTQ_ROOT_DIR = "dir";
    private static final String CSV_OUT_DIR = "out";
    private static final String THREAD_COUNT = "threadCount";

    public static void main(String[] args) throws ParseException, IOException, InterruptedException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String filePath = cmd.getOptionValue(FASTQ_FILE);
        final String dirPath = cmd.getOptionValue(FASTQ_ROOT_DIR);
        final String csvOutPath = cmd.getOptionValue(CSV_OUT_DIR);
        final String threadCountArg = cmd.getOptionValue(THREAD_COUNT);

        if ((filePath == null && dirPath == null) || csvOutPath == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Fastq-Stats", options);
        } else if (filePath != null) {
            final FastqTracker tracker = FastqStats.processFile(filePath);
            writeOutputToCSV("", tracker, csvOutPath);
        } else {
            final int threadCount = getThreadCount(threadCountArg);
            final String flowcellName = getFlowcellName(dirPath);
            final File baseCallsDir = getBaseCallsDir(dirPath);
            final long startTime = System.currentTimeMillis();
            final FastqTracker tracker = FastqStats.processDir(baseCallsDir, threadCount);
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
        options.addOption(THREAD_COUNT, true, "Number of max threads to use (only used when running on a directory).");
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
        for (final String sampleName : tracker.samples().keySet()) {
            for (final String laneName : tracker.samples().get(sampleName).keySet()) {
                final FastqData sampleLaneData = tracker.samples().get(sampleName).get(laneName);
                writer.write(sampleName + ", " + laneName + ", " + sampleLaneData.yield() + ", "
                        + sampleLaneData.q30() * 100.0 / sampleLaneData.yield() + "\n");
            }
        }
        writer.close();
        LOGGER.info("Written fastq qc data to " + csvOutPath);
    }

    /**
     * Extracts the flowcell name from the RunInfo.xml file which is assumed to be located in the dirPath directory.
     * The flowcell name is assumed to be a text node under the <Flowcell> element.
     *
     * @param dirPath path to the directory
     * @return The flowcell name, or "Unknown" if the flowcell name couldn't be extracted.
     */

    @NotNull
    @VisibleForTesting
    static String getFlowcellName(@NotNull final String dirPath) {
        final File runInfoXml = new File(dirPath + File.separator + "RunInfo.xml");
        final DocumentBuilderFactory domFactory = DocumentBuilderFactory.newInstance();
        try {
            final DocumentBuilder documentBuilder = domFactory.newDocumentBuilder();
            final Document xmlDoc = documentBuilder.parse(runInfoXml);
            final XPath xPath = XPathFactory.newInstance().newXPath();
            return (String) xPath.evaluate("//Flowcell/text()", xmlDoc, XPathConstants.STRING);
        } catch (Exception e) {
            LOGGER.warn("Could not extract flowcell name from " + runInfoXml.getPath());
            return "Unknown";
        }
    }

    @NotNull
    @VisibleForTesting
    static File getBaseCallsDir(@NotNull final String dirPath) throws IOException {
        final File baseCallsDir = new File(
                dirPath + File.separator + "Data" + File.separator + "Intensities" + File.separator + "BaseCalls");
        if (!baseCallsDir.exists()) {
            throw new IOException("dir " + baseCallsDir + " does not exist.");
        } else if (!baseCallsDir.isDirectory()) {
            throw new IOException(baseCallsDir + " is not a directory.");
        }
        return baseCallsDir;
    }

    @VisibleForTesting
    static int getThreadCount(@NotNull final String threadCountArg) {
        try {
            final int numThreads = Integer.parseInt(threadCountArg);
            if (numThreads <= 0) {
                throw new NumberFormatException();
            }
            return numThreads;
        } catch (NumberFormatException e) {
            LOGGER.info("Couldn't parse thread count parameter > 0; using default value.");
            return Runtime.getRuntime().availableProcessors();
        }
    }
}
