package com.hartwig.hmftools.pave.pon_gen;

import java.util.List;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

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

        int firstPosMatchIndex = -1;
        boolean matched = false;

        for(; mCurrentIndex < mHotspots.size(); ++mCurrentIndex)
        {
            VariantHotspot hotspot = mHotspots.get(mCurrentIndex);

            if(hotspot.position() < position)
                continue;

            if(hotspot.position() > position)
                break;

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
        else if(mCurrentIndex > 0)
            --mCurrentIndex;

        return matched;
    }
}

