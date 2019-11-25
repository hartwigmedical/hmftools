package com.hartwig.hmftools.sage.select;

import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.variant.SageVariantTier;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class TierSelector {

    private final RegionSelector<GenomeRegion> regionSelector;
    private final PositionSelector<VariantHotspot> hotspotSelector;

    public TierSelector(final List<GenomeRegion> panel, final List<VariantHotspot> hotspots) {
        this.regionSelector = new RegionSelector<>(panel);
        this.hotspotSelector = new PositionSelector<>(hotspots);
    }

    public boolean isHotspot(@NotNull final VariantHotspot variant) {
        final AtomicBoolean hotspotMatch = new AtomicBoolean(false);

        hotspotSelector.select(variant.position(), variant.position(), hotspot -> {
            if (hotspot.alt().equals(variant.alt()) && hotspot.ref().equals(variant.ref())) {
                hotspotMatch.set(true);
            }
        });

        return hotspotMatch.get();
    }

    @NotNull
    public SageVariantTier tier(@NotNull final VariantHotspot variant) {

        if (isHotspot(variant)) {
            return SageVariantTier.HOTSPOT;
        }

        if (regionSelector.select(variant.position()).isPresent()) {
            return SageVariantTier.PANEL;
        }

        return SageVariantTier.WIDE;
    }

}
