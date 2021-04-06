package com.hartwig.hmftools.cobalt.diploid;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CountBamLinesDiploidBed implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(CountBamLinesDiploidBed.class);

    private static final String DELIMITER = "\t";

    public static void main(String[] args) throws ParseException {
        final Options options = createOptions();
        final CommandLine cmd = new DefaultParser().parse(options, args);
        final String inputFile = cmd.getOptionValue("in");
        final String outputFile = cmd.getOptionValue("out");

        try (CountBamLinesDiploidBed app = new CountBamLinesDiploidBed(inputFile, outputFile)) {
            app.run();
        } catch (Exception e) {
            LOGGER.error(e);
            e.printStackTrace();
            System.exit(-1);
        }
    }

    private final String inputFile;
    private final String outputFile;
    private final long timestamp = System.currentTimeMillis();

    public CountBamLinesDiploidBed(final String inputFile, final String outputFile) {
        this.inputFile = inputFile;
        this.outputFile = outputFile;
    }

    public void run() throws IOException {
        double cutoff = 0.50;

        LOGGER.info("Reading input file: {}", inputFile);
        final List<DiploidCount> diploidCounts = DiploidCount.readDiploidCountAsList(inputFile);

        LOGGER.info("Determining diploid regions with {} cutoff", cutoff);
//        final DiploidRegionBuilder builder = new DiploidRegionBuilder(cutoff, 100, 100);
        final DiploidRegionBuilder builder = new DiploidRegionBuilder(cutoff, 62, 34);
        diploidCounts.forEach(builder);
        final List<GenomeRegion> diploidRegions = builder.build();
        LOGGER.info("Total diploid bases: {}", builder.getTotalDiploidBases());

        Files.write(new File(outputFile).toPath(),
                diploidRegions.stream().sorted().map(CountBamLinesDiploidBed::toBedString).collect(Collectors.toList()));
    }

    @NotNull
    private static String toBedString(@NotNull GenomeRegion region) {
        return new StringJoiner(DELIMITER).add(region.chromosome())
                .add(String.valueOf(region.start() - 1))
                .add(String.valueOf(region.end()))
                .toString();
    }

    @Override
    public void close() {
        LOGGER.info("Complete in {} seconds", (System.currentTimeMillis() - timestamp) / 1000);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption("in", true, "Cobalt diploid count");
        options.addOption("out", true, "Bed File");
        return options;
    }
}
