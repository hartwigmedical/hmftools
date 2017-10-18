package com.hartwig.hmftools.common.pcf;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

public class PCFFile {

    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "sampleID";
    private static final String RATIO_EXTENSION = ".cobalt.ratio.pcf";
    private static final String BAF_EXTENSION = ".amber.baf.pcf";

    @NotNull
    public static String generateRatioFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + RATIO_EXTENSION;
    }

    @NotNull
    public static String generateBAFFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + BAF_EXTENSION;
    }

    public static Multimap<Chromosome, PCFPosition> readPositions(@NotNull final String filename) throws IOException {
        Multimap<Chromosome, PCFPosition> result = ArrayListMultimap.create();
        for (String line : Files.readAllLines(new File(filename).toPath())) {
            if (!line.startsWith(HEADER_PREFIX)) {
                String[] values = line.split(DELIMITER);
                final String chromosomeName = values[1];
                if (HumanChromosome.contains(chromosomeName)) {
                    final Chromosome chromosome = HumanChromosome.fromString(chromosomeName);
                    result.put(chromosome, position(chromosomeName, Long.valueOf(values[3])));
                    result.put(chromosome, position(chromosomeName, Long.valueOf(values[4])));
                }
            }
        }

        return result;
    }

    private static PCFPosition position(@NotNull final String chromosome, long pos) {
        return  ImmutablePCFPosition.builder().chromosome(chromosome).position(pos).build();
    }

    @NotNull
    public static Multimap<String, PCFRegion> read(int windowSize, @NotNull final String filename) throws IOException {
        return fromLines(windowSize, Files.readAllLines(new File(filename).toPath()));
    }

    @NotNull
    private static Multimap<String, PCFRegion> fromLines(int windowSize, @NotNull List<String> lines) {
        Multimap<String, PCFRegion> result = ArrayListMultimap.create();
        for (String line : lines) {
            if (!line.startsWith(HEADER_PREFIX)) {
                final PCFRegion region = fromString(windowSize, line);
                result.put(region.chromosome(), region);
            }
        }
        return result;
    }

    @NotNull
    private static PCFRegion fromString(int windowSize, @NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return ImmutablePCFRegion.builder()
                .chromosome(values[1])
                .start(Long.valueOf(values[3]))
                .end(Long.valueOf(values[4]) + windowSize - 1)
                .build();
    }
}
