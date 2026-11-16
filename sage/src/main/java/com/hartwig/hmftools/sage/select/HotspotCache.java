package com.hartwig.hmftools.sage.select;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.variant.SimpleVariant;

public class HotspotCache
{
    // non-thread safe cache of hotspots
    private final List<SimpleVariant> mHotspots;
    private final int mMinPosition;
    private final int mMaxPosition;
    private int mIndex;

    public HotspotCache(final List<SimpleVariant> hotspots)
    {
        mHotspots = hotspots;
        mIndex = 0;

        if(!mHotspots.isEmpty())
        {
            mMinPosition = mHotspots.get(0).Position;
            mMaxPosition = mHotspots.get(mHotspots.size() - 1).Position;
        }
        else
        {
            mMaxPosition = -1;
            mMinPosition = -1;
        }
    }

    public boolean isHotspot(final SimpleVariant variant)
    {
        if(mHotspots.isEmpty())
            return false;

        if(!positionWithin(variant.Position, mMinPosition, mMaxPosition))
            return false;

        // move to first position match if already matching
        while(mIndex >= 0 && mHotspots.get(mIndex).Position == variant.Position)
        {
            --mIndex;

            if(mIndex < 0 || mHotspots.get(mIndex).Position < variant.Position)
            {
                ++mIndex;
                break;
            }
        }

        if(mHotspots.get(mIndex).Position <= variant.Position)
        {
            // move forward through the list
            if(mHotspots.get(mIndex).matches(variant))
                return true;

            ++mIndex;

            while(mIndex < mHotspots.size())
            {
                SimpleVariant hotspot = mHotspots.get(mIndex);

                if(hotspot.Position > variant.Position)
                    break;

                if(hotspot.matches(variant))
                    return true;

                ++mIndex;
            }
        }
        else
        {
            // move back through the list
            --mIndex;

            while(mIndex >= 0)
            {
                SimpleVariant hotspot = mHotspots.get(mIndex);

                if(hotspot.Position < variant.Position)
                    return false;

                if(hotspot.matches(variant))
                    return true;

                --mIndex;
            }
        }

        mIndex = min(max(0, mIndex), mHotspots.size() - 1);

        return false;
    }

    @VisibleForTesting
    public int index() { return mIndex; }
}
