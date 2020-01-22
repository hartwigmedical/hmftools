package com.hartwig.hmftools.sage.select;

import java.util.List;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.variant.SageVariantTier;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class TierSelector extends HotspotSelector{

    private final PanelSelector<GenomeRegion> regionSelector;

    public TierSelector(final List<GenomeRegion> panel, final List<VariantHotspot> hotspots) {
        super(hotspots);
        this.regionSelector = new PanelSelector<>(panel);
    }

    @NotNull
    public SageVariantTier tier(@NotNull final VariantHotspot variant) {

        if (isHotspot(variant)) {
            return SageVariantTier.HOTSPOT;
        }

        if (regionSelector.inPanel(variant.position(), variant.end())) {
            return SageVariantTier.PANEL;
        }

        return SageVariantTier.WIDE;
    }

}
