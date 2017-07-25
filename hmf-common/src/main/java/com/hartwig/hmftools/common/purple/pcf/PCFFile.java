package com.hartwig.hmftools.common.purple.pcf;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

import org.jetbrains.annotations.NotNull;

public class PCFFile {

    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "sampleID";
    private static final String EXTENSION = ".purple.pcf";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
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
                .start(Long.valueOf(values[3]) + windowSize - 1)
                .end(Long.valueOf(values[4]))
                .build();
    }
}
