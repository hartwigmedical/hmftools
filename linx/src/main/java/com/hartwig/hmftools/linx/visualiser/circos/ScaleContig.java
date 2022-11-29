package com.hartwig.hmftools.linx.visualiser.circos;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class ScaleContig
{
    private int mMaxAdjustedPosition;
    private final String mContig;
    private final Map<Integer,Integer> mPositionMap;

    @VisibleForTesting
    ScaleContig(@NotNull final String contig, Map<Integer, Integer> positionMap) 
    {
        mContig = contig;
        mPositionMap = positionMap;
    }

    ScaleContig(@NotNull final String contig, @NotNull final List<Integer> positions)
    {
        this.mContig = contig;
        this.mPositionMap = Maps.newHashMap();
        final List<Integer> sortedDistinctPositions = positions.stream().sorted().distinct().collect(Collectors.toList());

        if(!sortedDistinctPositions.isEmpty())
        {
            int currentAdjustedPosition = 1;
            int previousUnadjustedPositionPosition = sortedDistinctPositions.get(0);
            mPositionMap.put(previousUnadjustedPositionPosition, currentAdjustedPosition);

            for (int i = 1; i < sortedDistinctPositions.size(); i++)
            {
                int currentUnadjustedPosition = sortedDistinctPositions.get(i);
                int linearDistance = currentUnadjustedPosition - previousUnadjustedPositionPosition;
                int logDistance = logDistance(linearDistance);
                currentAdjustedPosition = currentAdjustedPosition + logDistance;
                mPositionMap.put(currentUnadjustedPosition, currentAdjustedPosition);

                previousUnadjustedPositionPosition = currentUnadjustedPosition;
            }

            mMaxAdjustedPosition = currentAdjustedPosition;
        }
    }

    public boolean isEmpty()
    {
        return mPositionMap.isEmpty();
    }

    public String contig()
    {
        return mContig;
    }

    public int length()
    {
        return mMaxAdjustedPosition;
    }

    public void expand(double factor)
    {
        for (Map.Entry<Integer, Integer> entry : mPositionMap.entrySet())
        {
            if(entry.getValue() > 1)
            {
                entry.setValue((int) Math.round(factor * entry.getValue()));
            }
        }

        mMaxAdjustedPosition = (int) Math.round(factor * mMaxAdjustedPosition);
    }

    public int scale(int position)
    {
        if(!mPositionMap.containsKey(position))
        {
            throw new IllegalArgumentException("Invalid position");
        }

        return mPositionMap.get(position);
    }

    public int interpolate(int value)
    {
        final Set<Integer> keySet = mPositionMap.keySet();

        if(mPositionMap.containsKey(value))
        {
            return mPositionMap.get(value);
        }

        int minValue = keySet.stream().mapToInt(x -> x).min().orElse(0);
        int maxValue = keySet.stream().mapToInt(x -> x).max().orElse(0);

        int closestToStart = keySet.stream().filter(x -> x < value).mapToInt(x -> x).max().orElse(minValue);
        int closestToEnd = keySet.stream().filter(x -> x > value).mapToInt(x -> x).min().orElse(maxValue);
        if(closestToStart == closestToEnd)
        {
            return mPositionMap.get(closestToStart);
        }

        double longDistanceProportion = Math.abs(value - closestToStart) / ((double) Math.abs(closestToEnd - closestToStart));

        int closestIntToStart = mPositionMap.get(closestToStart);
        int closestIntToEnd = mPositionMap.get(closestToEnd);

        return closestIntToStart + (int) Math.floor(longDistanceProportion * Math.abs(closestIntToEnd - closestIntToStart));
    }

    static int logDistance(int distance)
    {
        return (int) Math.floor(Math.pow(Math.log10(distance), 3)) + 10;
    }

}
