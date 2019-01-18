package com.hartwig.hmftools.svanalysis.visualisation;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;

import org.jetbrains.annotations.NotNull;

public class ScalePosition {

    private final Map<String, Map<Long, Integer>> chromosomePositionMap = Maps.newHashMap();

    public ScalePosition(final int start, @NotNull final List<GenomeRegion> regions) {
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
    public List<GenomeRegion> scaleRegions(@NotNull final List<GenomeRegion> regions) {
        return regions.stream().map(x -> scale(x, chromosomePositionMap.get(x.chromosome()))).collect(Collectors.toList());
    }

    @NotNull
    public List<SvLink> scaleLinks(@NotNull final List<SvLink> links) {
        return links.stream()
                .map(x -> scale(x, chromosomePositionMap.get(x.startChromosome()), chromosomePositionMap.get(x.endChromosome())))
                .collect(Collectors.toList());
    }

    @NotNull
    private static GenomeRegion scale(@NotNull final GenomeRegion victim, @NotNull final Map<Long, Integer> positionMap) {
        return GenomeRegionFactory.create(victim.chromosome(), positionMap.get(victim.start()), positionMap.get(victim.end()));
    }

    @NotNull
    private static SvLink scale(@NotNull final SvLink victim, @NotNull final Map<Long, Integer> startPositionMap,
            @NotNull final Map<Long, Integer> endPositionMap) {
        return new SvLink(victim.startChromosome(),
                startPositionMap.get(victim.startPosition()),
                victim.endChromosome(),
                endPositionMap.get(victim.endPosition()));
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
