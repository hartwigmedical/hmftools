package com.hartwig.hmftools.sage.select;

import java.util.List;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.variant.SageVariantTier;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class TierSelector extends HotspotSelector
{
    private final PanelSelector<GenomeRegion> mPanelRegionSelector;
    private final PanelSelector<GenomeRegion> mHighConfidenceRegionSelector;

    public TierSelector(final List<VariantHotspot> hotspots, final List<GenomeRegion> panel, final List<GenomeRegion> highConfidence)
    {
        super(hotspots);
        this.mPanelRegionSelector = new PanelSelector<>(panel);
        this.mHighConfidenceRegionSelector = new PanelSelector<>(highConfidence);
    }

    @NotNull
    public SageVariantTier tier(@NotNull final VariantHotspot variant)
    {
        if(isHotspot(variant))
        {
            return SageVariantTier.HOTSPOT;
        }

        if(mPanelRegionSelector.inPanel(variant.position(), variant.end()))
        {
            return SageVariantTier.PANEL;
        }

        if(mHighConfidenceRegionSelector.inPanel(variant.position(), variant.end()))
        {
            return SageVariantTier.HIGH_CONFIDENCE;
        }

        return SageVariantTier.LOW_CONFIDENCE;
    }
}
