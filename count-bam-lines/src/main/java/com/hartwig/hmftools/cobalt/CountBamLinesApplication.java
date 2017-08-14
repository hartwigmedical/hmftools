package com.hartwig.hmftools.cobalt;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;

import com.google.common.collect.Lists;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.chromosome.ChromosomeLengthFactory;
import com.hartwig.hmftools.common.chromosome.ChromosomeLengthFile;
import com.hartwig.hmftools.common.cobalt.ReadCount;
import com.hartwig.hmftools.common.cobalt.ReadCountFile;

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

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class CountBamLinesApplication {
    private static final Logger LOGGER = LogManager.getLogger(CountBamLinesApplication.class);

    private static final String THREADS = "threads";
    private static final String INPUT_FILE = "input";
    private static final String OUTPUT_DIR = "output_dir";
    private static final String WINDOW_SIZE = "window_size";
    private static final String SAMPLE = "sample";
    private static final String MIN_QUALITY = "min_quality";

    private static final int WINDOW_SIZE_DEFAULT = 1000;
    private static final int MIN_QUALITY_DEFAULT = 10;

    public static void main(final String... args) throws ParseException, IOException, ExecutionException, InterruptedException {
        new CountBamLinesApplication(args);
    }

    private CountBamLinesApplication(final String... args) throws ParseException, IOException, ExecutionException, InterruptedException {

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        if (!cmd.hasOption(INPUT_FILE) || !cmd.hasOption(OUTPUT_DIR) || !cmd.hasOption(SAMPLE)) {
            printUsageAndExit(options);
        }

        final String sample = cmd.getOptionValue(SAMPLE);
        final Path outputPath = selectOrCreateOutputPath(cmd.getOptionValue(OUTPUT_DIR));
        if (outputPath == null) {
            System.exit(1);
        }

        final File inputFile = new File(cmd.getOptionValue(INPUT_FILE));
        final String outputCountFile = ReadCountFile.generateFilename(outputPath.toString(), sample);
        final String chromosomeLengthFile = ChromosomeLengthFile.generateFilename(outputPath.toString(), sample);
        final int threadCount = cmd.hasOption(THREADS) ? Integer.valueOf(cmd.getOptionValue(THREADS)) : 4;
        final int windowSize = cmd.hasOption(WINDOW_SIZE) ? Integer.valueOf(cmd.getOptionValue(WINDOW_SIZE)) : WINDOW_SIZE_DEFAULT;
        final int minQuality = cmd.hasOption(MIN_QUALITY) ? Integer.valueOf(cmd.getOptionValue(MIN_QUALITY)) : MIN_QUALITY_DEFAULT;
        LOGGER.info("Input File: {}", inputFile.toString());
        LOGGER.info("Output Read Counts: {}", outputCountFile);
        LOGGER.info("Output Chromosome Lengths: {}", chromosomeLengthFile);
        LOGGER.info("Thread Count: {}, Window Size: {}, Min Quality {}", threadCount, windowSize, minQuality);

        final SamReaderFactory readerFactory = SamReaderFactory.make();
        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("bam-%d").build();
        final ExecutorService executorService = Executors.newFixedThreadPool(threadCount, namedThreadFactory);

        final List<ChromosomeLength> lengths;
        try (SamReader reader = readerFactory.open(inputFile)) {
            lengths = ChromosomeLengthFactory.create(reader.getFileHeader());
        }
        ChromosomeLengthFile.write(chromosomeLengthFile, lengths);

        final List<Future<ChromosomeReadCount>> futures = Lists.newArrayList();
        for (ChromosomeLength chromosome : lengths) {
            final ChromosomeReadCount callable = new ChromosomeReadCount(inputFile,
                    readerFactory,
                    chromosome.chromosome(),
                    chromosome.position(),
                    windowSize,
                    minQuality);
            futures.add(executorService.submit(callable));
        }

        ReadCountFile.createFile(windowSize, outputCountFile);
        for (Future<ChromosomeReadCount> future : futures) {
            persist(outputCountFile, future.get());
        }

        LOGGER.info("Complete");
        executorService.shutdown();
    }

    @Nullable
    private Path selectOrCreateOutputPath(@NotNull final String outputString) {
        final File outputFile = new File(outputString);
        if (!outputFile.exists() && !outputFile.mkdirs()) {
            LOGGER.error("Unable to create output directory {} ", outputString);
            return null;
        }

        if (!outputFile.isDirectory()) {
            LOGGER.error("Output dir {} is not a valid directory", outputString);
            return null;
        }

        return outputFile.toPath();
    }

    private void persist(@NotNull final String filename, @NotNull ChromosomeReadCount counter) throws IOException {
        final List<ReadCount> readCounts = counter.readCount();
        LOGGER.info("Persisting {} windows from chromosome {}", readCounts.size(), counter.chromosome());
        ReadCountFile.append(filename, readCounts);
    }

    private static void printUsageAndExit(@NotNull final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("COBALT", options);
        System.exit(1);
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        try {
            return parser.parse(options, args);
        } catch (ParseException e) {
            printUsageAndExit(options);
            throw e;
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(WINDOW_SIZE, true, "Window size. Default 1000.");
        options.addOption(THREADS, true, "Number of threads. Default 4.");
        options.addOption(INPUT_FILE, true, "Input bam location/filename");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(SAMPLE, true, "Name of sample");

        return options;
    }
}
