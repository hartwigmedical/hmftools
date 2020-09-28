package com.hartwig.hmftools.serve.hotspot;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HotspotFunctions {

    private static final Logger LOGGER = LogManager.getLogger(HotspotFunctions.class);

    private HotspotFunctions() {
    }

    @NotNull
    public static Map<VariantHotspot, HotspotAnnotation> convertHotspotMap(@NotNull String source,
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
                    LOGGER.debug("Annotation '{}' already found previously: '{}'", newAnnotation, currentAnnotation);
                    convertedMap.put(hotspot, mergeHotspotAnnotations(currentAnnotation, newAnnotation));
                } else {
                    convertedMap.put(hotspot, newAnnotation);
                }
            }
        }

        return convertedMap;
    }

    @NotNull
    public static Map<VariantHotspot, HotspotAnnotation> mergeHotspotMaps(@NotNull List<Map<VariantHotspot, HotspotAnnotation>> maps) {
        Map<VariantHotspot, HotspotAnnotation> mergedMap = Maps.newTreeMap(new VariantHotspotComparator());
        for (Map<VariantHotspot, HotspotAnnotation> map : maps) {
            for (Map.Entry<VariantHotspot, HotspotAnnotation> entry : map.entrySet()) {
                VariantHotspot hotspot = entry.getKey();
                HotspotAnnotation newAnnotation = entry.getValue();
                HotspotAnnotation currentAnnotation = mergedMap.get(hotspot);

                if (currentAnnotation != null) {
                    mergedMap.put(hotspot, mergeHotspotAnnotations(currentAnnotation, newAnnotation));
                } else {
                    mergedMap.put(hotspot, newAnnotation);
                }
            }
        }
        return mergedMap;
    }

    @NotNull
    public static HotspotAnnotation mergeHotspotAnnotations(@NotNull HotspotAnnotation annotation1,
            @NotNull HotspotAnnotation annotation2) {
        if (!annotation1.gene().equals(annotation2.gene())) {
            LOGGER.warn("Genes mismatch on identical hotspot: '{}' vs '{}'", annotation1.gene(), annotation2.gene());
        }

        String bestTranscript = annotation1.transcript();
        String bestProteinAnnotation = annotation1.proteinAnnotation();
        if (bestTranscript == null) {
            bestTranscript = annotation2.transcript();
            bestProteinAnnotation = annotation2.proteinAnnotation();
        }

        Set<String> mergedSources = Sets.newHashSet();
        mergedSources.addAll(annotation1.sources());
        mergedSources.addAll(annotation2.sources());

        return new HotspotAnnotation(mergedSources, annotation1.gene(), bestTranscript, bestProteinAnnotation);
    }
}
