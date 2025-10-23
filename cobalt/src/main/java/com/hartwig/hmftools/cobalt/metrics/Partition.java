package com.hartwig.hmftools.cobalt.metrics;

import static java.lang.String.format;

import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class Partition extends ChrBaseRegion
{
    public final List<TargetRegionData> TargetRegions;
    private final int WindowSize;

    public Partition(String chromosome, int start, int end, int windowSize)
    {
        super(chromosome, start, end);
        int length = end - start + 1;
        Preconditions.checkArgument(length >= windowSize);
        Preconditions.checkArgument(end % 1000 == 0);
        Preconditions.checkArgument(length % 1000 == 0);
        Preconditions.checkArgument(windowSize % 1000 == 0);
        int windowCount = length / windowSize;
        Preconditions.checkArgument(windowCount % 10 == 0);
        TargetRegions = new ArrayList<>(windowCount);
        WindowSize = windowSize;
        for(int i = 0; i < windowCount; i++)
        {
            int s = i * windowSize;
            TargetRegions.add(new TargetRegionData(chromosome, start + s, start + s + windowSize - 1));
        }
    }

    public void recordFragment(int startPosition, int length)
    {
        if(!containsPosition(startPosition))
        {
            // e.g. negative read with alignment start near partition start.
            return;
        }
        Preconditions.checkArgument(containsPosition(startPosition));
        int offset = startPosition - start();
        int index = offset / WindowSize;
        TargetRegions.get(index).recordFragment(length);
    }

    public String toString()
    {
        return format("region(%s) targeted(%d)", super.toString(), TargetRegions.size());
    }
}
