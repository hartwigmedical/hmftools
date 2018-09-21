package com.hartwig.hmftools.common.variant.enrich;

import java.util.Collection;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class HotspotEnrichment implements SomaticEnrichment {

    private static final int DISTANCE = 5;
    private static final String HOTSPOT_TAG = "HOTSPOT";

    private final Multimap<String, VariantHotspot> hotspots;

    public HotspotEnrichment(@NotNull final Multimap<String, VariantHotspot> hotspots) {
        this.hotspots = hotspots;
    }

    @NotNull
    @Override
    public ImmutableSomaticVariantImpl.Builder enrich(@NotNull final ImmutableSomaticVariantImpl.Builder builder,
            @NotNull final VariantContext context) {
        if (context.hasAttribute(HOTSPOT_TAG)) {
            builder.hotspot(Hotspot.HOTSPOT);
        } else if (hotspots.containsKey(context.getContig()) && overlaps(context, hotspots.get(context.getContig()))) {
            builder.hotspot(Hotspot.NEAR_HOTSPOT);
        } else {
            builder.hotspot(Hotspot.NON_HOTSPOT);
        }

        return builder;
    }

    private static boolean overlaps(@NotNull final VariantContext variant, Collection<VariantHotspot> hotspots) {
        return hotspots.stream().anyMatch(x -> overlaps(x, variant));
    }

    @VisibleForTesting
    static boolean overlaps(@NotNull final VariantHotspot hotspot, @NotNull final VariantContext variant) {
        int variantStart = variant.getStart();
        int variantEnd = variant.getStart() + variant.getReference().length() - 1 + DISTANCE;

        long ponStart = hotspot.position();
        long ponEnd = hotspot.position() + hotspot.ref().length() - 1 + DISTANCE;

        return variantStart <= ponEnd && variantEnd >= ponStart;
    }
}
