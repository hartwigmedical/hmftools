package com.hartwig.hmftools.svanalysis.visualisation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;

import org.jetbrains.annotations.NotNull;

public class Links {

    private static final String HEADER = "SampleId";
    private static final String COMMENT = "#";
    private static final String DELIMITER = ",";

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

    public static int linkTraverseCount(@NotNull final GenomePosition position, @NotNull final List<Link> links) {
       return  Links.findStartLink(position, links).map(Link::traverseCount).orElse(0) + Links.findEndLink(position, links)
                .map(Link::traverseCount)
                .orElse(0);
    }

    @NotNull
    public static List<Link> clean(@NotNull final List<Link> links) {
        return links.stream().filter(x -> HumanChromosome.contains(x.startChromosome()) && HumanChromosome.contains(x.endChromosome())).collect(Collectors.toList());
    }

    @NotNull
    public static List<Link> readLinks(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    private static List<Link> fromLines(@NotNull List<String> lines) {
        final List<Link> results = Lists.newArrayList();
        for (final String line : lines) {
            if (!line.startsWith(COMMENT) && !line.startsWith(HEADER) ) {
                results.add(fromString(line));
            }
        }

        return results;
    }

    @NotNull
    private static Link fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableLink.builder()
                .clusterId(Integer.valueOf(values[1]))
                .chainId(Integer.valueOf(values[2]))
                .startChromosome(values[4])
                .startPosition(Long.valueOf(values[5]))
                .startOrientation(Integer.valueOf(values[6]))
                .startType(Link.Type.valueOf(values[7]))
                .endChromosome(values[8])
                .endPosition(Long.valueOf(values[9]))
                .endOrientation(Integer.valueOf(values[10]))
                .endType(Link.Type.valueOf(values[11]))
                .traverseCount(Integer.valueOf(values[12]))
                .build();
    }

    @NotNull
    public static List<GenomePosition> allPositions(@NotNull final List<Link> links) {
        final List<GenomePosition> results = Lists.newArrayList();

        for (final Link link : links) {
            if (HumanChromosome.contains(link.startChromosome()) && link.startPosition() != -1) {
                results.add(GenomePositions.create(link.startChromosome(), link.startPosition()));
            }

            if (HumanChromosome.contains(link.endChromosome()) && link.endPosition() != -1) {
                results.add(GenomePositions.create(link.endChromosome(), link.endPosition()));
            }
        }

        Collections.sort(results);

        return results;
    }

}
