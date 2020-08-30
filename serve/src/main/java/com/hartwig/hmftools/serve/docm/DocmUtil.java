package com.hartwig.hmftools.serve.docm;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.serve.hotspot.HotspotAnnotation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class DocmUtil {

    private static final Logger LOGGER = LogManager.getLogger(DocmUtil.class);

    private DocmUtil() {
    }

    @NotNull
    public static Map<VariantHotspot, HotspotAnnotation> convertHotspots(@NotNull Map<DocmEntry, List<VariantHotspot>> hotspotsPerEntry) {
        Map<VariantHotspot, HotspotAnnotation> convertedMap = Maps.newTreeMap(new VariantHotspotComparator());
        for (Map.Entry<DocmEntry, List<VariantHotspot>> entry : hotspotsPerEntry.entrySet()) {
            DocmEntry docmEntry = entry.getKey();
            for (VariantHotspot hotspot : entry.getValue()) {
                HotspotAnnotation newAnnotation = new HotspotAnnotation(Sets.newHashSet("docm"),
                        docmEntry.gene(),
                        docmEntry.transcript(),
                        docmEntry.proteinAnnotation());

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
