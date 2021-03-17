package com.hartwig.hmftools.serve.extraction.hotspot;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotImpl;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class HotspotFunctions {

    private static final Logger LOGGER = LogManager.getLogger(HotspotFunctions.class);

    private HotspotFunctions() {
    }

    @NotNull
    public static Set<KnownHotspot> consolidate(@NotNull Iterable<KnownHotspot> hotspots) {
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
                LOGGER.debug("Merging hotspots {} with {} on {}", newAnnotation, existingAnnotation, key);
                annotationPerHotspot.put(key, mergeHotspotAnnotations(newAnnotation, existingAnnotation));
            }
        }

        Set<KnownHotspot> consolidatedHotspots = Sets.newHashSet();
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

        String bestTranscript;
        String bestProteinAnnotation;

        // If both annotations either have or have no transcript annotation it does not matter which one we pick, but we do enforce a
        // choice to make sure hotspot files stay identical for identical inputs.
        boolean favorAnnotation1;
        if (annotation1.proteinAnnotation().equals(annotation2.proteinAnnotation())) {
            if (annotation1.transcript() != null && annotation2.transcript() != null) {
                favorAnnotation1 = annotation1.transcript().compareTo(annotation2.transcript()) > 0;
            } else {
                favorAnnotation1 = annotation1.transcript() != null;
            }
        } else {
            favorAnnotation1 = annotation1.proteinAnnotation().compareTo(annotation2.proteinAnnotation()) > 0;
        }

        if (annotation1.transcript() == null && annotation2.transcript() == null) {
            bestTranscript = null;
            bestProteinAnnotation = favorAnnotation1 ? annotation1.proteinAnnotation() : annotation2.proteinAnnotation();
        } else if (annotation1.transcript() == null) {
            bestTranscript = annotation2.transcript();
            bestProteinAnnotation = annotation2.proteinAnnotation();
        } else if (annotation2.transcript() == null) {
            bestTranscript = annotation1.transcript();
            bestProteinAnnotation = annotation1.proteinAnnotation();
        } else {
            bestTranscript = favorAnnotation1 ? annotation1.transcript() : annotation2.transcript();
            bestProteinAnnotation = favorAnnotation1 ? annotation1.proteinAnnotation() : annotation2.proteinAnnotation();
        }

        Set<Knowledgebase> mergedSources = Sets.newHashSet();
        mergedSources.addAll(annotation1.sources());
        mergedSources.addAll(annotation2.sources());

        return new HotspotAnnotation(mergedSources, annotation1.gene(), bestTranscript, bestProteinAnnotation);
    }

    private static class HotspotAnnotation {

        @NotNull
        private final Set<Knowledgebase> sources;
        @NotNull
        private final String gene;
        @Nullable
        private final String transcript;
        @NotNull
        private final String proteinAnnotation;

        public HotspotAnnotation(@NotNull final Set<Knowledgebase> sources, @NotNull final String gene, @Nullable final String transcript,
                @NotNull final String proteinAnnotation) {
            this.sources = sources;
            this.gene = gene;
            this.transcript = transcript;
            this.proteinAnnotation = proteinAnnotation;
        }

        @NotNull
        public Set<Knowledgebase> sources() {
            return sources;
        }

        @NotNull
        public String gene() {
            return gene;
        }

        @Nullable
        public String transcript() {
            return transcript;
        }

        @NotNull
        public String proteinAnnotation() {
            return proteinAnnotation;
        }

        @Override
        public String toString() {
            return "HotspotAnnotation{" + "sources=" + sources + ", gene='" + gene + '\'' + ", transcript='" + transcript + '\''
                    + ", proteinAnnotation='" + proteinAnnotation + '\'' + '}';
        }
    }
}
