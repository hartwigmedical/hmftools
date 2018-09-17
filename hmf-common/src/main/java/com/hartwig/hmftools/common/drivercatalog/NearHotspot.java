package com.hartwig.hmftools.common.drivercatalog;

import java.util.Collection;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class NearHotspot {

    private final Multimap<String, GenomePosition> hotspots;
    private final int distance;

    public NearHotspot(int distance, @NotNull final Multimap<String, GenomePosition> hotspots) {
        this.distance = distance;
        this.hotspots = hotspots;
    }

    public NearHotspot(int distance, List<GenomePosition> hotspots) {
        this.distance = distance;
        this.hotspots = ArrayListMultimap.create();
        for (GenomePosition hotspot : hotspots) {
            this.hotspots.put(hotspot.chromosome(), hotspot);
        }
    }

    boolean isNearHotspot(@NotNull final SomaticVariant variant) {
        return hotspots.containsKey(variant.chromosome()) && isNearHotspot(variant, hotspots.get(variant.chromosome()));
    }

    private boolean isNearHotspot(@NotNull final SomaticVariant variant, Collection<GenomePosition> hotspots) {
        long start = variant.position() - distance;
        long end = variant.position() + variant.ref().length() - 1 + distance;

        for (GenomePosition hotspot : hotspots) {
            if (hotspot.position() >= start && hotspot.position() <= end) {
                return true;
            }
        }

        return false;
    }

}
