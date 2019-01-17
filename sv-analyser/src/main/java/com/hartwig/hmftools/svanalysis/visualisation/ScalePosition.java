package com.hartwig.hmftools.svanalysis.visualisation;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class ScalePosition {

    private static final int BUFFER = 1;

    private long minPosition;

    public ScalePosition(long minPosition) {
        this.minPosition = minPosition;
    }

    @NotNull
    public List<GenomeRegion> scale(@NotNull final List<GenomeRegion> regions) {
        return regions;
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
