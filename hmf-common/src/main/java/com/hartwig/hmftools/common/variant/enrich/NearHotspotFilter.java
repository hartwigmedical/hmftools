package com.hartwig.hmftools.common.variant.enrich;

import java.util.Collection;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class NearHotspotFilter implements VariantContextFilter {

    private final Multimap<String, VariantHotspot> hotspots;
    private final int distance;

    public NearHotspotFilter(int distance, @NotNull final Multimap<String, VariantHotspot> hotspots) {
        this.distance = distance;
        this.hotspots = hotspots;
    }

    public NearHotspotFilter(int distance, @NotNull final List<VariantHotspot> hotspots) {
        this.distance = distance;
        this.hotspots = ArrayListMultimap.create();
        for (VariantHotspot hotspot : hotspots) {
            this.hotspots.put(hotspot.chromosome(), hotspot);
        }
    }

    @Override
    public boolean test(final VariantContext variant) {
        return hotspots.containsKey(variant.getContig()) && isNearHotspot(variant, hotspots.get(variant.getContig()));
    }

    private boolean isNearHotspot(@NotNull final VariantContext variant, Collection<VariantHotspot> hotspots) {
        long start = variant.getStart();
        long end = variant.getStart() + variant.getReference().getBaseString().length() - 1 + distance;

        for (VariantHotspot hotspot : hotspots) {
            if (hotspot.position() >= start && hotspot.position() + hotspot.ref().length() - 1 + distance <= end) {
                return true;
            }
        }

        return false;
    }

}
