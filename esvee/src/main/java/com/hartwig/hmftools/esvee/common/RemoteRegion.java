package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class RemoteRegion
{
    private ChrBaseRegion mRegion;
    private final Set<String> mReadIds;

    public RemoteRegion(final ChrBaseRegion region, final String readId)
    {
        mRegion = region;
        mReadIds = Sets.newHashSet(readId);
    }

    public void addReadDetails(final String readId, final int posStart, final int posEnd)
    {
        mRegion.setStart(min(mRegion.start(), posStart));
        mRegion.setEnd(max(mRegion.end(), posEnd));
        mReadIds.add(readId);
    }

    public boolean overlaps(final String chromosome, final int posStart, final int posEnd)
    {
        return mRegion.overlaps(chromosome, posStart, posEnd);
    }

    public ChrBaseRegion region() { return mRegion; }

    public Set<String> readIds() { return mReadIds; }
    public int readCount() { return mReadIds.size(); }

    public String toString() { return format("%s reads(%d)", mRegion, mReadIds.size()); }

    public static void mergeRegions(final List<RemoteRegion> regions)
    {
        Collections.sort(regions, Comparator.comparing(x -> x.region()));

        int index = 0;
        while(index < regions.size() - 1)
        {
            RemoteRegion region = regions.get(index);

            int nextIndex = index + 1;
            while(nextIndex < regions.size())
            {
                RemoteRegion nextRegion = regions.get(nextIndex);

                if(region.region().overlaps(nextRegion.region()))
                {
                    regions.remove(nextIndex);
                    region.region().setEnd(nextRegion.region().end());
                    region.readIds().addAll(nextRegion.readIds());
                }
                else
                {
                    ++nextIndex;
                }
            }

            ++index;
        }
    }
}
