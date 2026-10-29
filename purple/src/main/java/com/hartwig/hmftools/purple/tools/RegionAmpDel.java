package com.hartwig.hmftools.purple.tools;

import java.util.Objects;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency;

import org.jetbrains.annotations.NotNull;

public class RegionAmpDel implements Comparable<RegionAmpDel>
{
    final ChrBaseRegion mRegion;
    final AmpDelRegionFrequency.EventType mType;

    public RegionAmpDel(RegionGeneEvents event)
    {
        mRegion = event.region();
        mType = event.eventType();
    }

    public RegionAmpDel(ChrBaseRegion region, AmpDelRegionFrequency.EventType type)
    {
        mRegion = region;
        mType = type;
    }

    public int start()
    {
        return mRegion.start();
    }

    public int end()
    {
        return mRegion.end();
    }

    public String type()
    {
        return mType.toString();
    }

    @Override
    public int compareTo(@NotNull final RegionAmpDel o)
    {
        int result = new EventRegionComparator().compare(mRegion, o.mRegion);
        if(result != 0)
        {
            return result;
        }
        result = mType.compareTo(o.mType);
        return result;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final RegionAmpDel that = (RegionAmpDel) o;
        return Objects.equals(mRegion, that.mRegion) && mType == that.mType;
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mRegion, mType);
    }

    @Override
    public String toString()
    {
        return String.format("%s:%d-%d:%s", mRegion.chromosome(), mRegion.start(), mRegion.end(), mType);
    }
}
