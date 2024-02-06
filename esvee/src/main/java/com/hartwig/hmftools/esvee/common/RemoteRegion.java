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

public class RemoteRegion extends ChrBaseRegion
{
    private final Set<String> mReadIds;

    public static final int REMOTE_READ_TYPE_JUNCTION_MATE = 0;
    public static final int REMOTE_READ_TYPE_JUNCTION_SUPP = 1;
    public static final int REMOTE_READ_TYPE_DISCORDANT_READ = 2;

    private final int[] mReadTypeCount;

    public RemoteRegion(final ChrBaseRegion region, final String readId, final int readType)
    {
        super(region.Chromosome, region.start(), region.end());
        mReadIds = Sets.newHashSet(readId);
        mReadTypeCount = new int[REMOTE_READ_TYPE_DISCORDANT_READ+1];
        ++mReadTypeCount[readType];
    }

    public void addReadDetails(final String readId, final int posStart, final int posEnd, final int readType)
    {
        setStart(min(start(), posStart));
        setEnd(max(end(), posEnd));
        mReadIds.add(readId);
        ++mReadTypeCount[readType];
    }

    public Set<String> readIds() { return mReadIds; }
    public int readCount() { return mReadIds.size(); }

    public int[] readTypeCounts() { return mReadTypeCount; }

    public String toString() { return format("%s reads(%d)", super.toString(), mReadIds.size()); }

    public static void mergeRegions(final List<RemoteRegion> regions)
    {
        Collections.sort(regions, Comparator.comparing(x -> x));

        int index = 0;
        while(index < regions.size() - 1)
        {
            RemoteRegion region = regions.get(index);

            int nextIndex = index + 1;
            while(nextIndex < regions.size())
            {
                RemoteRegion nextRegion = regions.get(nextIndex);

                if(region.overlaps(nextRegion))
                {
                    regions.remove(nextIndex);
                    region.setEnd(max(region.end(), nextRegion.end()));
                    region.readIds().addAll(nextRegion.readIds());

                    for(int i = 0; i < region.readTypeCounts().length; ++i)
                    {
                        region.readTypeCounts()[i] += nextRegion.readTypeCounts()[i];
                    }
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
