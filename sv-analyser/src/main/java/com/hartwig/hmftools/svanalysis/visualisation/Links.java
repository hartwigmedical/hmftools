package com.hartwig.hmftools.svanalysis.visualisation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

public class Links {

    private static final String COMMENT = "#";
    private static final String DELIMITER = "\t";


    @NotNull
    static Optional<Link> findStartLink(@NotNull final GenomePosition position, @NotNull List<Link> links) {
        return links.stream()
                .filter(x -> x.startChromosome().equals(position.chromosome()) && x.startPosition() == position.position())
                .findFirst();
    }

    @NotNull
    static Optional<Link> findEndLink(@NotNull final GenomePosition position, @NotNull List<Link> links) {
        return links.stream()
                .filter(x -> x.endChromosome().equals(position.chromosome()) && x.endPosition() == position.position())
                .findFirst();
    }

    @NotNull
    public static Optional<Link> findLink(@NotNull final GenomePosition position, @NotNull final List<Link> links) {
        final Optional<Link> result = findStartLink(position, links);
        return result.isPresent() ? result : findEndLink(position, links);
    }

    @NotNull
    public static List<Link> clean(@NotNull final List<Link> links) {
        return links.stream().filter(x -> x.startPosition() != -1 && x.endPosition() != -1).collect(Collectors.toList());
    }

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
