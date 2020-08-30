package com.hartwig.hmftools.serve.hotspot;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HotspotUtil {

    private static final Logger LOGGER = LogManager.getLogger(HotspotUtil.class);

    private HotspotUtil() {
    }

    @NotNull
    public static Map<VariantHotspot, HotspotAnnotation> convertHotspots(@NotNull String source,
            @NotNull Map<? extends HotspotSourceEntry, List<VariantHotspot>> hotspotsPerEntry) {
        Map<VariantHotspot, HotspotAnnotation> convertedMap = Maps.newTreeMap(new VariantHotspotComparator());
        for (Map.Entry<? extends HotspotSourceEntry, List<VariantHotspot>> entry : hotspotsPerEntry.entrySet()) {
            HotspotSourceEntry sourceEntry = entry.getKey();
            for (VariantHotspot hotspot : entry.getValue()) {
                HotspotAnnotation newAnnotation = new HotspotAnnotation(Sets.newHashSet(source),
                        sourceEntry.gene(),
                        sourceEntry.transcript(),
                        sourceEntry.proteinAnnotation());

                HotspotAnnotation currentAnnotation = convertedMap.get(hotspot);
                if (currentAnnotation != null) {
                    LOGGER.warn("Annotation '{}' already found previously: '{}'", newAnnotation, currentAnnotation);
                } else {
                    convertedMap.put(hotspot, newAnnotation);
                }
            }
        }

        return convertedMap;
    }
}
