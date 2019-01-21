package com.hartwig.hmftools.svanalysis.visualisation;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class ScalePosition {

    private final Map<String, Map<Long, Integer>> chromosomePositionMap = Maps.newHashMap();

    public ScalePosition(final int start, @NotNull final List<? extends GenomeRegion> regions) {
        final Set<String> contigs = regions.stream().map(GenomeRegion::chromosome).collect(Collectors.toSet());
        for (final String contig : contigs) {
            final List<Long> contigPositions = Lists.newArrayList();
            regions.stream().filter(x -> x.chromosome().equals(contig)).forEach(x -> {
                contigPositions.add(x.start());
                contigPositions.add(x.end());
            });

            chromosomePositionMap.put(contig, positionMap(start, contigPositions));
        }
    }

    @NotNull
    public List<Track> scaleTracks(@NotNull final List<Track> tracks) {
        return tracks.stream().map(x -> scale(x, chromosomePositionMap.get(x.chromosome()))).collect(Collectors.toList());
    }

    @NotNull
    public List<Link> scaleLinks(@NotNull final List<Link> links) {
        return links.stream()
                .map(x -> scale(x, chromosomePositionMap.get(x.startChromosome()), chromosomePositionMap.get(x.endChromosome())))
                .collect(Collectors.toList());
    }

    @NotNull
    private static Track scale(@NotNull final Track victim, @NotNull final Map<Long, Integer> positionMap) {
        return ImmutableTrack.builder()
                .chromosome(victim.chromosome())
                .start(positionMap.get(victim.start()))
                .end(positionMap.get(victim.end()))
                .track(victim.track())
                .build();
    }

    @NotNull
    private static Link scale(@NotNull final Link victim, @NotNull final Map<Long, Integer> startPositionMap,
            @NotNull final Map<Long, Integer> endPositionMap) {

        return ImmutableLink.builder()
                .startChromosome(victim.startChromosome())
                .startPosition(startPositionMap.get(victim.startPosition()))
                .endChromosome(victim.endChromosome())
                .endPosition(endPositionMap.get(victim.endPosition()))
                .build();

    }

    @VisibleForTesting
    static Map<Long, Integer> positionMap(int start, @NotNull final Long... positionArray) {
        return positionMap(start, Lists.newArrayList(positionArray));
    }

    @NotNull
    private static Map<Long, Integer> positionMap(int start, @NotNull final List<Long> positions) {
        final Map<Long, Integer> results = Maps.newHashMap();
        final List<Long> sortedDistinctPositions = positions.stream().sorted().distinct().collect(Collectors.toList());

        if (!sortedDistinctPositions.isEmpty()) {
            int logPosition = start;
            long lastPosition = sortedDistinctPositions.get(0);
            results.put(lastPosition, logPosition);

            for (int i = 1; i < sortedDistinctPositions.size(); i++) {
                long position = sortedDistinctPositions.get(i);
                long linearDistance = position - lastPosition;
                int logDistance = logDistance(linearDistance);
                logPosition = logPosition + logDistance;
                lastPosition = position;

                results.put(lastPosition, logPosition);
            }
        }

        return results;
    }

    static int logDistance(long distance) {
        return (int) Math.floor(Math.pow(Math.log10(distance), 2)) + 1;
    }

}
