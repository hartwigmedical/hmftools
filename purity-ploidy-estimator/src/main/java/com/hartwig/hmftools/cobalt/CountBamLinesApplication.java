package com.hartwig.hmftools.cobalt;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;

import com.google.common.collect.Lists;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.chromosome.ChromosomeLengthFile;
import com.hartwig.hmftools.common.cobalt.ReadCount;
import com.hartwig.hmftools.common.cobalt.ReadCountFile;
import com.hartwig.hmftools.purple.LoadSomaticVariants;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReaderFactory;

public class CountBamLinesApplication {
    private static final Logger LOGGER = LogManager.getLogger(LoadSomaticVariants.class);

    private static final String THREADS = "threads";
    private static final String BAM_FILE = "bam_file";
    private static final String COUNT_FILE = "count_file";
    private static final String CHR_LENGTHS = "chr_len_file";
    private static final int WINDOW_SIZE = 1000;

    public static void main(final String... args) throws ParseException, IOException, ExecutionException, InterruptedException {
        new CountBamLinesApplication(args);
    }

    private CountBamLinesApplication(final String... args) throws ParseException, IOException, ExecutionException, InterruptedException {

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        if (!cmd.hasOption(BAM_FILE) || !cmd.hasOption(COUNT_FILE) || !cmd.hasOption(CHR_LENGTHS)) {
            printUsageAndExit(options);
        }

        final File inputFile = new File(cmd.getOptionValue(BAM_FILE));
        final String outputFile = cmd.getOptionValue(COUNT_FILE);
        final int threadCount = cmd.hasOption(THREADS) ? Integer.valueOf(cmd.getOptionValue(THREADS)) : 4;

        final SamReaderFactory readerFactory = SamReaderFactory.make();
        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("bam-%d").build();
        final ExecutorService executorService = Executors.newFixedThreadPool(threadCount, namedThreadFactory);

        final List<Future<ChromosomeCount>> futures = Lists.newArrayList();
        final List<ChromosomeLength> lengths = ChromosomeLengthFile.read(cmd.getOptionValue(CHR_LENGTHS));
        for (ChromosomeLength chromosome : lengths) {
            final CallableChromosome callable =
                    new CallableChromosome(chromosome.chromosome(), chromosome.position(), WINDOW_SIZE, readerFactory, inputFile);
            futures.add(executorService.submit(callable));
        }

        ReadCountFile.createFile(WINDOW_SIZE, outputFile);
        for (Future<ChromosomeCount> future : futures) {
            persist(outputFile, future.get());
        }

        LOGGER.info("Complete");
        executorService.shutdown();
    }

    private void persist(@NotNull final String filename, @NotNull ChromosomeCount counter) throws IOException {
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
        return parser.parse(options, args);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(THREADS, true, "Number of threads. Default 4.");
        options.addOption(CHR_LENGTHS, true, "Location of chromosome lengths file");
        options.addOption(BAM_FILE, true, "Location of input bam file");
        options.addOption(COUNT_FILE, true, "Location of output count file");

        return options;
    }
}
