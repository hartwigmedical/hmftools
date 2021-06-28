package com.hartwig.hmftools.sage.select;

import java.util.List;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.variant.VariantTier;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class TierSelector extends HotspotSelector
{
    private final PanelSelector mPanelRegionSelector;
    private final PanelSelector mHighConfidenceRegionSelector;

    public TierSelector(final List<VariantHotspot> hotspots, final List<BaseRegion> panel, final List<BaseRegion> highConfidence)
    {
        super(hotspots);
        mPanelRegionSelector = new PanelSelector(panel);
        mHighConfidenceRegionSelector = new PanelSelector(highConfidence);
    }

    @NotNull
    public VariantTier tier(@NotNull final VariantHotspot variant)
    {
        if(isHotspot(variant))
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
