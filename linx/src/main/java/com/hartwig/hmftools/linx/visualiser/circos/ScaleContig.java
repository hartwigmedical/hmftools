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
    private final Map<Integer, Integer> positionMap;


    @VisibleForTesting
    ScaleContig(@NotNull final String contig, Map<Integer, Integer> positionMap) {
        this.contig = contig;
        this.positionMap = positionMap;
    }

    ScaleContig(@NotNull final String contig, @NotNull final List<Integer> positions)
    {
        this.contig = contig;
        this.positionMap = Maps.newHashMap();
        final List<Integer> sortedDistinctPositions = positions.stream().sorted().distinct().collect(Collectors.toList());

        if (!sortedDistinctPositions.isEmpty())
        {
            int currentAdjustedPosition = 1;
            int previousUnadjustedPositionPosition = sortedDistinctPositions.get(0);
            positionMap.put(previousUnadjustedPositionPosition, currentAdjustedPosition);

            for (int i = 1; i < sortedDistinctPositions.size(); i++)
            {
                int currentUnadjustedPosition = sortedDistinctPositions.get(i);
                int linearDistance = currentUnadjustedPosition - previousUnadjustedPositionPosition;
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
        for (Map.Entry<Integer, Integer> entry : positionMap.entrySet())
        {
            if (entry.getValue() > 1)
            {
                entry.setValue((int) Math.round(factor * entry.getValue()));
            }
        }

        maxAdjustedPosition = (int) Math.round(factor * maxAdjustedPosition);
    }

    public int scale(int position)
    {
        if (!positionMap.containsKey(position))
        {
            throw new IllegalArgumentException("Invalid position");
        }

        return positionMap.get(position);
    }

    public int interpolate(int value)
    {

        final Set<Integer> keySet = positionMap.keySet();

        if (positionMap.containsKey(value))
        {
            return positionMap.get(value);
        }

        int minValue = keySet.stream().mapToInt(x -> x).min().orElse(0);
        int maxValue = keySet.stream().mapToInt(x -> x).max().orElse(0);

        int closestToStart = keySet.stream().filter(x -> x < value).mapToInt(x -> x).max().orElse(minValue);
        int closestToEnd = keySet.stream().filter(x -> x > value).mapToInt(x -> x).min().orElse(maxValue);
        if (closestToStart == closestToEnd)
        {
            return positionMap.get(closestToStart);
        }

        double longDistanceProportion = Math.abs(value - closestToStart) / ((double) Math.abs(closestToEnd - closestToStart));

        int closestIntToStart = positionMap.get(closestToStart);
        int closestIntToEnd = positionMap.get(closestToEnd);

        return closestIntToStart + (int) Math.floor(longDistanceProportion * Math.abs(closestIntToEnd - closestIntToStart));
    }

    static int logDistance(int distance)
    {
        return (int) Math.floor(Math.pow(Math.log10(distance), 3)) + 10;
    }

}
