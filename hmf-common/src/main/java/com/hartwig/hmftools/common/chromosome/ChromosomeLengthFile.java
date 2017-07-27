package com.hartwig.hmftools.common.chromosome;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class ChromosomeLengthFile {

    private static final String DELIMITER = "\t";

    @NotNull
    public static Map<String, ChromosomeLength> read(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    @NotNull
    private static Map<String, ChromosomeLength> fromLines(@NotNull List<String> lines) {
        Map<String, ChromosomeLength> result = Maps.newHashMap();
        for (String line : lines) {
            final ChromosomeLength length = fromString(line);
            result.put(length.chromosome(), length);
        }
        return result;
    }

    @NotNull
    private static ChromosomeLength fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableChromosomeLength.builder().chromosome(values[0]).position(Long.valueOf(values[2])).build();
    }

}
