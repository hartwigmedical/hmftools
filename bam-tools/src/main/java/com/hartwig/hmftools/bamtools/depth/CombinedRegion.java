package com.hartwig.hmftools.bamtools.depth;

import static java.lang.Math.max;
import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;

class CombinedRegion
{
    public List<PositionCount> Depth;

    public CombinedRegion(final HighDepthRegion region)
    {
        Depth = Lists.newArrayList();

        for(int pos = region.start(); pos <= region.end(); ++pos)
        {
            Depth.add(new PositionCount(pos, region.DepthMin, region.DepthMax));
        }
    }

    public int start()
    {
        return Depth.get(0).Position;
    }

    public int end()
    {
        return Depth.get(Depth.size() - 1).Position;
    }

    public int length()
    {
        return Depth.size();
    }

    public void addBases(final HighDepthRegion region)
    {
        for(int pos = region.start(); pos <= region.end(); ++pos)
        {
            int existingIndex = 0;

            boolean found = false;
            while(existingIndex < Depth.size())
            {
                PositionCount existing = Depth.get(existingIndex);

                if(pos == existing.Position)
                {
                    found = true;
                    ++existing.Count;
                    existing.DepthMax = max(existing.DepthMax, region.DepthMax);
                    break;
                }

                if(pos < existing.Position)
                {
                    break;
                }

                ++existingIndex;
            }

            if(!found)
            {
                Depth.add(existingIndex, new PositionCount(pos, region.DepthMin, region.DepthMax));
            }
        }
    }

    public void addRegion(final CombinedRegion other)
    {
        for(PositionCount otherCount : other.Depth)
        {
            int existingIndex = 0;

            boolean found = false;
            while(existingIndex < Depth.size())
            {
                PositionCount existing = Depth.get(existingIndex);

                if(otherCount.Position == existing.Position)
                {
                    found = true;
                    existing.Count += otherCount.Count;
                    break;
                }

                if(otherCount.Position < existing.Position)
                {
                    break;
                }

                ++existingIndex;
            }

            if(!found)
            {
                Depth.add(existingIndex, otherCount);
            }
        }
    }

    public String toString()
    {
        return format("span(%d - %d) length(%d)", start(), end(), length());
    }
}
