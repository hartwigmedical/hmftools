package com.hartwig.hmftools.linx.visualiser.circos;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

class ScaleContig
{
    private int maxAdjustedPosition;
    private final String contig;
    private final Map<Long, Integer> positionMap;


    @VisibleForTesting
    ScaleContig(@NotNull final String contig, Map<Long, Integer> positionMap) {
        this.contig = contig;
        this.positionMap = positionMap;
    }

    ScaleContig(@NotNull final String contig, @NotNull final List<Long> positions)
    {
        this.contig = contig;
        this.positionMap = Maps.newHashMap();
        final List<Long> sortedDistinctPositions = positions.stream().sorted().distinct().collect(Collectors.toList());

        if (!sortedDistinctPositions.isEmpty())
        {
            int currentAdjustedPosition = 1;
            long previousUnadjustedPositionPosition = sortedDistinctPositions.get(0);
            positionMap.put(previousUnadjustedPositionPosition, currentAdjustedPosition);

            for (int i = 1; i < sortedDistinctPositions.size(); i++)
            {
                long currentUnadjustedPosition = sortedDistinctPositions.get(i);
                long linearDistance = currentUnadjustedPosition - previousUnadjustedPositionPosition;
                int logDistance = logDistance(linearDistance);
                currentAdjustedPosition = currentAdjustedPosition + logDistance;
                positionMap.put(currentUnadjustedPosition, currentAdjustedPosition);

                previousUnadjustedPositionPosition = currentUnadjustedPosition;
            }

            maxAdjustedPosition = currentAdjustedPosition;
        }

    }

    public boolean isEmpty()
    {
        return positionMap.isEmpty();
    }

    public String contig()
    {
        return contig;
    }

    public int length()
    {
        return maxAdjustedPosition;
    }

    public void expand(double factor)
    {
        for (Map.Entry<Long, Integer> entry : positionMap.entrySet())
        {
            if (entry.getValue() > 1)
            {
                entry.setValue((int) Math.round(factor * entry.getValue()));
            }
        }

        maxAdjustedPosition = (int) Math.round(factor * maxAdjustedPosition);
    }

    public int scale(long position)
    {
        if (!positionMap.containsKey(position))
        {
            throw new IllegalArgumentException("Invalid position");
        }

        return positionMap.get(position);
    }

    public int interpolate(long value)
    {

        final Set<Long> keySet = positionMap.keySet();

        if (positionMap.containsKey(value))
        {
            return positionMap.get(value);
        }

        long minValue = keySet.stream().mapToLong(x -> x).min().orElse(0);
        long maxValue = keySet.stream().mapToLong(x -> x).max().orElse(0);

        long closestToStart = keySet.stream().filter(x -> x < value).mapToLong(x -> x).max().orElse(minValue);
        long closestToEnd = keySet.stream().filter(x -> x > value).mapToLong(x -> x).min().orElse(maxValue);
        if (closestToStart == closestToEnd)
        {
            return positionMap.get(closestToStart);
        }

        double longDistanceProportion = Math.abs(value - closestToStart) / ((double) Math.abs(closestToEnd - closestToStart));

        int closestIntToStart = positionMap.get(closestToStart);
        int closestIntToEnd = positionMap.get(closestToEnd);

        return closestIntToStart + (int) Math.floor(longDistanceProportion * Math.abs(closestIntToEnd - closestIntToStart));
    }

    static int logDistance(long distance)
    {
        return (int) Math.floor(Math.pow(Math.log10(distance), 3)) + 10;
    }

}
