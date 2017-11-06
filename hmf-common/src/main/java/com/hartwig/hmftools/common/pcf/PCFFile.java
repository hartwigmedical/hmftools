package com.hartwig.hmftools.common.pcf;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.window.Window;

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

    public static ListMultimap<String, PCFPosition> readPositions(int windowSize, @NotNull PCFSource source, @NotNull final String filename)
            throws IOException {
        ListMultimap<String, PCFPosition> result = ArrayListMultimap.create();
        final Window window = new Window(windowSize);
        long end = 0;
        for (String line : Files.readAllLines(new File(filename).toPath())) {
            if (!line.startsWith(HEADER_PREFIX)) {
                String[] values = line.split(DELIMITER);
                final String chromosomeName = values[1];
                long start = window.start(Long.valueOf(values[3]));
                if (start != end) {
                    result.put(chromosomeName, position(chromosomeName, start, source));
                }
                end = window.start(Long.valueOf(values[4])) + windowSize;
                result.put(chromosomeName, position(chromosomeName, end, source));
            }
        }

        return result;
    }

    @NotNull
    public static Multimap<String, GenomeRegion> read(int windowSize, @NotNull final String filename) throws IOException {
        return fromLines(windowSize, Files.readAllLines(new File(filename).toPath()));
    }

    @NotNull
    private static Multimap<String, GenomeRegion> fromLines(int windowSize, @NotNull List<String> lines) {
        Multimap<String, GenomeRegion> result = ArrayListMultimap.create();
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

    private static PCFPosition position(@NotNull final String chromosome, long pos, @NotNull final PCFSource source) {
        return ImmutablePCFPosition.builder().chromosome(chromosome).position(pos).source(source).build();
    }

}
