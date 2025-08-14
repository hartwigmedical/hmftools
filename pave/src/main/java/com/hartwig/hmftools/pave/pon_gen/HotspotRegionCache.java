package com.hartwig.hmftools.pave.pon_gen;

import static java.lang.String.format;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.variant.VariantHotspot;

public class HotspotRegionCache
{
    private final List<VariantHotspot> mHotspots;
    private int mCurrentIndex;

    public HotspotRegionCache(final List<VariantHotspot> hotspots)
    {
        mHotspots = hotspots;
        mCurrentIndex = 0;
    }

    public void resetSearch() { mCurrentIndex = 0; }
    public int entryCount() { return mHotspots.size(); }

    public boolean matchesHotspot(final int position, final String ref, final String alt)
    {
        if(mHotspots == null || mHotspots.isEmpty())
            return false;

        // cannot search at earlier positions than the current index
        if(mCurrentIndex < mHotspots.size() && position < mHotspots.get(mCurrentIndex).position())
            return false;

        int firstPosMatchIndex = -1;
        boolean matched = false;

        for(; mCurrentIndex < mHotspots.size(); ++mCurrentIndex)
        {
            VariantHotspot hotspot = mHotspots.get(mCurrentIndex);

            if(hotspot.position() < position)
                continue;

            if(hotspot.position() > position)
            {
                if(firstPosMatchIndex == -1)
                {
                    // no match on position was found - move index back to prior hotspot
                    if(mCurrentIndex > 0)
                        --mCurrentIndex;
                }
                else
                {
                    // otherwise set to first match at this position ready for the next call / check
                    mCurrentIndex = firstPosMatchIndex;
                }

                return false;
            }

            if(firstPosMatchIndex == -1)
                firstPosMatchIndex = mCurrentIndex;

            if(hotspot.ref().equals(ref) && hotspot.alt().equals(alt))
            {
                matched = true;
                break;
            }
        }

        // move the index back to the prior position or the first at this position
        if(firstPosMatchIndex >= 0)
            mCurrentIndex = firstPosMatchIndex;

        return matched;
    }

    @VisibleForTesting
    public int currentIndex() { return mCurrentIndex; }

    public String toString()
    {
        return format("size(%d) index(%d) currentPos(%d)",
            mHotspots.size(), mCurrentIndex, mCurrentIndex < mHotspots.size() ? mHotspots.get(mCurrentIndex).position() : -1);
    }
}

