package com.hartwig.hmftools.svanalysis.visualisation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class Links {

    private static final String COMMENT = "#";
    private static final String DELIMITER = "\t";

    @NotNull
    public static List<Link> readLinks(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    private static List<Link> fromLines(@NotNull List<String> lines) {
        final List<Link> results = Lists.newArrayList();
        for (final String line : lines) {
            if (!line.startsWith(COMMENT)) {
                results.add(fromString(line));
            }
        }

        return results;
    }

    @NotNull
    private static Link fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableLink.builder()
                .chainId(Integer.valueOf(values[0]))
                .startChromosome(values[1])
                .startPosition(Long.valueOf(values[2]))
                .startOrientation(Integer.valueOf(values[3]))
                .startFoldback(Boolean.valueOf(values[4]))
                .endChromosome(values[5])
                .endPosition(Long.valueOf(values[6]))
                .endOrientation(Integer.valueOf(values[7]))
                .endFoldback(Boolean.valueOf(values[8]))
                .build();
    }

}
