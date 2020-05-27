package com.hartwig.hmftools.serve.vicc.range;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.jetbrains.annotations.NotNull;

public class GeneLevelEventExtractor {

    private static final Set<String> ONCOKB_GENE_ACTIVATION =
            Sets.newHashSet("Gain-of-function Mutations");

    private static final Set<String> ONCOKB_GENE_INACTIVATION =
            Sets.newHashSet("Truncating Mutations");

    @NotNull
    public Map<Feature, String> extractKnownGeneLevelEvents(@NotNull ViccEntry viccEntry) {
        Map<Feature, String> geneLevelEventsPerFeature = Maps.newHashMap();
        if (viccEntry.source() == ViccSource.ONCOKB) {
            for (Feature feature : viccEntry.features()) {
                if (ONCOKB_GENE_ACTIVATION.contains(feature.name())) {
                    geneLevelEventsPerFeature.put(feature, "gain of " + feature.geneSymbol());
                } else if (ONCOKB_GENE_INACTIVATION.contains(feature.name())) {
                    geneLevelEventsPerFeature.put(feature, "loss of " + feature.geneSymbol());
                }
            }
        }
        return geneLevelEventsPerFeature;
    }
}
