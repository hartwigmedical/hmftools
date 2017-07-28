package com.hartwig.hmftools.cobalt;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.chromosome.ChromosomeLengthFile;
import com.hartwig.hmftools.common.purple.ratio.ReadCount;
import com.hartwig.hmftools.common.purple.ratio.ReadCountFile;
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
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class CountBamLinesApplication {
    private static final Logger LOGGER = LogManager.getLogger(LoadSomaticVariants.class);

    private static final String CHR_LENGTHS = "chr_len_file";
    private static final String BAM_FILE = "bam_file";
    private static final String COUNT_FILE = "count_file";
    private static final int WINDOW_SIZE = 1000;

    public static void main(final String... args) throws ParseException, IOException {
        new CountBamLinesApplication(args);
    }

    private CountBamLinesApplication(final String... args) throws ParseException, IOException {

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        if (!cmd.hasOption(BAM_FILE) || !cmd.hasOption(COUNT_FILE)) {
            printUsageAndExit(options);
        }

        final Map<String, ChromosomeLength> chromosomeLength = chromosomeLength(cmd);

        final File inputFile = new File(cmd.getOptionValue(BAM_FILE));
        final String outputFile = cmd.getOptionValue(COUNT_FILE);
        ReadCountFile.createFile(outputFile);

        SamReaderFactory readerFactory = SamReaderFactory.make();
        SamReader reader = readerFactory.open(inputFile);
        SAMRecordIterator iterator = reader.iterator();

        ChromosomeCount counter = null;
        while (iterator.hasNext()) {

            SAMRecord record = iterator.next();
            String contig = record.getContig();
            if (contig != null) {

                if (counter == null || !counter.contig().equals(contig)) {
                    persist(outputFile, counter);

                    final String chromosome = chromosome(contig);
                    final long length = Optional.ofNullable(chromosomeLength.get(chromosome)).map(ChromosomeLength::position).orElse(0L);

                    LOGGER.info("Generating windows on chromosome {}", chromosome);
                    counter = new ChromosomeCount(contig, chromosome, length, WINDOW_SIZE);
                }
                counter.addRecord(record);
            }
        }

        persist(outputFile, counter);
        LOGGER.info("Complete");
    }

    private void persist(@NotNull final String filename, @Nullable ChromosomeCount counter) throws IOException {
        if (counter != null) {
            final List<ReadCount> readCounts = counter.readCount();
            LOGGER.info("Persisting {} windows from chromosome {}", readCounts.size(), counter.chromosome());
            ReadCountFile.append(filename, readCounts);
        }
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
        options.addOption(CHR_LENGTHS, true, "Location of chromosome lengths file");
        options.addOption(BAM_FILE, true, "Location of input bam file");
        options.addOption(COUNT_FILE, true, "Location of output count file");

        return options;
    }

    @NotNull
    private Map<String, ChromosomeLength> chromosomeLength(@NotNull final CommandLine cmd) throws IOException {
        return cmd.hasOption(CHR_LENGTHS) ? ChromosomeLengthFile.read(cmd.getOptionValue(CHR_LENGTHS)) : Maps.newHashMap();
    }

    @NotNull
    private static String chromosome(@NotNull final String contig) {
        final String lower = contig.toLowerCase();
        return lower.startsWith("chr") ? contig.substring(3) : contig;
    }
}
