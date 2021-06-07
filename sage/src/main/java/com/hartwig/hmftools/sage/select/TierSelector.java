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

    private final PanelSelector<GenomeRegion> panelRegionSelector;
    private final PanelSelector<GenomeRegion> highConfidenceRegionSelector;

    public TierSelector(final List<VariantHotspot> hotspots, final List<GenomeRegion> panel, final List<GenomeRegion> highConfidence)
    {
        super(hotspots);
        this.panelRegionSelector = new PanelSelector<>(panel);
        this.highConfidenceRegionSelector = new PanelSelector<>(highConfidence);
    }

    @NotNull
    public SageVariantTier tier(@NotNull final VariantHotspot variant)
    {
        if(isHotspot(variant))
        {
            return SageVariantTier.HOTSPOT;
        }

        if(panelRegionSelector.inPanel(variant.position(), variant.end()))
        {
            return SageVariantTier.PANEL;
        }

        if(highConfidenceRegionSelector.inPanel(variant.position(), variant.end()))
        {
            return SageVariantTier.HIGH_CONFIDENCE;
        }

        return SageVariantTier.LOW_CONFIDENCE;
    }
}
