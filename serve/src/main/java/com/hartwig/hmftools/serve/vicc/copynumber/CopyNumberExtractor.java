package com.hartwig.hmftools.serve.vicc.copynumber;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.annotation.FeatureType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;

public class CopyNumberExtractor {

    @NotNull
    private final Set<String> uniqueAmps = Sets.newHashSet();
    @NotNull
    private final Set<String> uniqueDels = Sets.newHashSet();

    @NotNull
    public Set<String> uniqueAmps() {
        return uniqueAmps;
    }

    @NotNull
    public Set<String> uniqueDels() {
        return uniqueDels;
    }

    @NotNull
    public Map<Feature, KnownAmplificationDeletion> extractKnownAmplificationsDeletions(@NotNull ViccEntry viccEntry) {
        Map<Feature, KnownAmplificationDeletion> ampsDelsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            if (feature.type() == FeatureType.AMPLIFICATION) {
                ampsDelsPerFeature.put(feature, eventForGene(feature.geneSymbol(), CopyNumberType.AMPLIFICATION));
                uniqueAmps.add(feature.geneSymbol());
            } else if (feature.type() == FeatureType.DELETION) {
                ampsDelsPerFeature.put(feature, eventForGene(feature.geneSymbol(), CopyNumberType.DELETION));
                uniqueDels.add(feature.geneSymbol());
            }
        }
        return ampsDelsPerFeature;
    }

    @NotNull
    private static KnownAmplificationDeletion eventForGene(@NotNull String gene, @NotNull CopyNumberType type) {
        return ImmutableKnownAmplificationDeletion.builder().gene(gene).type(type).build();
    }
}
