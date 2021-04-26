package com.hartwig.hmftools.fastqstats;

import static com.hartwig.hmftools.fastqstats.FastqStats.getFastqsFromBaseCallsDir;
import static com.hartwig.hmftools.fastqstats.FastqStats.getFastqsFromDir;
import static com.hartwig.hmftools.fastqstats.FastqStats.getSingleFastq;

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
import com.google.common.collect.Multimap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.w3c.dom.Document;

public final class FastqStatsRunner {

    private static final Logger LOGGER = LogManager.getLogger(FastqStatsRunner.class);

    private static final String FASTQ_FILE = "file";
    private static final String FLOWCELL_ROOT_DIR = "dir";
    private static final String FASTQ_DIR = "fastq_dir";
    private static final String CSV_OUT_DIR = "out";
    private static final String THREAD_COUNT = "threadCount";

    public static void main(String[] args) throws ParseException, IOException, InterruptedException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String filePath = cmd.getOptionValue(FASTQ_FILE);
        final String flowcellDirPath = cmd.getOptionValue(FLOWCELL_ROOT_DIR);
        final String csvOutPath = cmd.getOptionValue(CSV_OUT_DIR);
        final String threadCountArg = cmd.getOptionValue(THREAD_COUNT);
        final String fastqDirPath = cmd.getOptionValue(FASTQ_DIR);

        if ((filePath == null && flowcellDirPath == null && fastqDirPath == null) || csvOutPath == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Fastq-Stats", options);
        } else if (filePath != null) {
            final Multimap<String, File> fastqsPerSample = getSingleFastq(filePath);
            final FastqTracker tracker = FastqStats.processFastqs(fastqsPerSample, 1);
            writeOutputToCSV("", tracker, csvOutPath);
        } else if (flowcellDirPath != null) {
            final int threadCount = getThreadCount(threadCountArg);
            final String flowcellName = getFlowcellName(flowcellDirPath);
            final File baseCallsDir = getBaseCallsDir(flowcellDirPath);
            final long startTime = System.currentTimeMillis();
            final Multimap<String, File> fastqsPerSample = getFastqsFromBaseCallsDir(baseCallsDir);
            final FastqTracker tracker = FastqStats.processFastqs(fastqsPerSample, threadCount);
            LOGGER.info("Total time: {}ms.", System.currentTimeMillis() - startTime);
            writeOutputToCSV(flowcellName, tracker, csvOutPath);
        } else {
            final int threadCount = getThreadCount(threadCountArg);
            final File fastqDir = getDir(fastqDirPath);
            final long startTime = System.currentTimeMillis();
            final Multimap<String, File> fastqsPerSample = getFastqsFromDir(fastqDir);
            final FastqTracker tracker = FastqStats.processFastqs(fastqsPerSample, threadCount);
            LOGGER.info("Total time: {}ms.", System.currentTimeMillis() - startTime);
            writeFastqDirOutputToCSV(tracker, csvOutPath);
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(FASTQ_FILE, true, "Path towards the original fastq file.");
        options.addOption(FLOWCELL_ROOT_DIR, true, "Path towards the flowcell dir.");
        options.addOption(CSV_OUT_DIR, true, "Path towards the csv output file.");
        options.addOption(THREAD_COUNT, true, "Number of max threads to use (only used when running on a directory).");
        options.addOption(FASTQ_DIR, true, "Path towards the fastq dir.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static void writeOutputToCSV(@NotNull final String flowcellName, @NotNull final FastqTracker tracker,
            @NotNull final String csvOutPath) throws IOException {
        final BufferedWriter writer = new BufferedWriter(new FileWriter(csvOutPath, false));
        writer.write("Flowcell," + flowcellName + "," + tracker.flowcell().yield() + "," + tracker.flowcell().q30Percentage() + "\n");

        for (final String laneName : tracker.lanes().keySet()) {
            final FastqData lane = tracker.lane(laneName);
            writer.write("Lane," + laneName + "," + lane.yield() + "," + lane.q30Percentage() + "\n");
        }
        for (final String sampleName : tracker.samples().keySet()) {
            final FastqData sample = tracker.sample(sampleName);
            writer.write("Sample," + sampleName + "," + sample.yield() + "," + sample.q30Percentage() + "\n");
        }
        for (final String sampleName : tracker.samples().keySet()) {
            for (final String laneName : tracker.samples().get(sampleName).keySet()) {
                final FastqData sampleLaneData = tracker.samples().get(sampleName).get(laneName);
                writer.write(sampleName + "," + laneName + "," + sampleLaneData.yield() + "," + sampleLaneData.q30Percentage() + "\n");
            }
        }
        writer.close();
        LOGGER.info("Written fastq qc data to {}", csvOutPath);
    }

    private static void writeFastqDirOutputToCSV(@NotNull final FastqTracker tracker, @NotNull final String csvOutPath) throws IOException {
        final BufferedWriter writer = new BufferedWriter(new FileWriter(csvOutPath, false));
        for (final String sampleName : tracker.samples().keySet()) {
            for (final String laneName : tracker.samples().get(sampleName).keySet()) {
                final FastqData sampleLaneData = tracker.samples().get(sampleName).get(laneName);
                writer.write(sampleName + "," + laneName + "," + sampleLaneData.yield() + "," + sampleLaneData.q30Percentage() + "\n");
            }
        }
        writer.close();
        LOGGER.info("Written fastq qc data to {}", csvOutPath);
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
            LOGGER.warn("Could not extract flowcell name from {}", runInfoXml.getPath());
            return "Unknown";
        }
    }

    @NotNull
    @VisibleForTesting
    static File getBaseCallsDir(@NotNull final String dirPath) throws IOException {
        return getDir(dirPath + File.separator + "Data" + File.separator + "Intensities" + File.separator + "BaseCalls");
    }

    @NotNull
    @VisibleForTesting
    static File getDir(@NotNull final String dirPath) throws IOException {
        final File dir = new File(dirPath);
        if (!dir.exists()) {
            throw new IOException("dir " + dirPath + " does not exist.");
        } else if (!dir.isDirectory()) {
            throw new IOException(dirPath + " is not a directory.");
        }
        return dir;
    }

    @VisibleForTesting
    static int getThreadCount(@Nullable final String threadCountArg) {
        final int availableThreads = Runtime.getRuntime().availableProcessors();
        if (threadCountArg == null) {
            return availableThreads;
        }

        try {
            final int numThreads = Integer.parseInt(threadCountArg);
            if (numThreads <= 0) {
                throw new NumberFormatException();
            }
            return numThreads;
        } catch (NumberFormatException e) {
            LOGGER.info("Couldn't parse thread count parameter > 0; using default value: {}.", availableThreads);
            return availableThreads;
        }
    }
}
