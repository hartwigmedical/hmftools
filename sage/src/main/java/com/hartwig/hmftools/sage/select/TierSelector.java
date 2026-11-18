package com.hartwig.hmftools.sage.select;

import java.util.List;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.common.variant.VariantTier;

@NotThreadSafe
public class TierSelector
{
    private final HotspotCache mHotspotCache;
    private final PanelSelector mPanelRegionSelector;
    private final PanelSelector mHighConfidenceRegionSelector;

    public TierSelector(final List<SimpleVariant> hotspots, final List<BaseRegion> panel, final List<BaseRegion> highConfidence)
    {
        mHotspotCache = new HotspotCache(hotspots);
        mPanelRegionSelector = new PanelSelector(panel);
        mHighConfidenceRegionSelector = new PanelSelector(highConfidence);
    }

    public VariantTier tier(final SimpleVariant variant)
    {
        if(mHotspotCache.isHotspot(variant))
        {
            return VariantTier.HOTSPOT;
        }

        if(mPanelRegionSelector.inPanel(variant.position(), variant.end()))
        {
            return VariantTier.PANEL;
        }

        if(mHighConfidenceRegionSelector.inPanel(variant.position(), variant.end()))
        {
            return VariantTier.HIGH_CONFIDENCE;
        }

        return VariantTier.LOW_CONFIDENCE;
    }
}
