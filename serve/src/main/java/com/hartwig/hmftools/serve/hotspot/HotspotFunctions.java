package com.hartwig.hmftools.serve.hotspot;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotImpl;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HotspotFunctions {

    private static final Logger LOGGER = LogManager.getLogger(HotspotFunctions.class);

    private HotspotFunctions() {
    }

    @NotNull
    public static List<KnownHotspot> consolidateHotspots(@NotNull List<KnownHotspot> hotspots) {
        Map<VariantHotspot, HotspotAnnotation> annotationPerHotspot = Maps.newHashMap();
        for (KnownHotspot hotspot : hotspots) {
            HotspotAnnotation newAnnotation = new HotspotAnnotation(Sets.newHashSet(hotspot.sources()),
                    hotspot.gene(),
                    hotspot.transcript(),
                    hotspot.proteinAnnotation());
            VariantHotspotImpl key = ImmutableVariantHotspotImpl.builder().from(hotspot).build();
            HotspotAnnotation existingAnnotation = annotationPerHotspot.get(key);
            if (existingAnnotation == null) {
                annotationPerHotspot.put(key, newAnnotation);
            } else {
                annotationPerHotspot.put(key, mergeHotspotAnnotations(newAnnotation, existingAnnotation));
            }
        }

        List<KnownHotspot> consolidatedHotspots = Lists.newArrayList();
        for (Map.Entry<VariantHotspot, HotspotAnnotation> entry : annotationPerHotspot.entrySet()) {
            HotspotAnnotation annotation = entry.getValue();
            consolidatedHotspots.add(ImmutableKnownHotspot.builder()
                    .from(entry.getKey())
                    .sources(annotation.sources())
                    .gene(annotation.gene())
                    .transcript(annotation.transcript())
                    .proteinAnnotation(annotation.proteinAnnotation())
                    .build());
        }
        return consolidatedHotspots;
    }

    @NotNull
    private static HotspotAnnotation mergeHotspotAnnotations(@NotNull HotspotAnnotation annotation1,
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
