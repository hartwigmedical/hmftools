package com.hartwig.hmftools.serve.vicc.signatures;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.jetbrains.annotations.NotNull;

public class SignaturesExtractor {

    private static final Set<String> ONCOKB_SIGNATURES =
            Sets.newHashSet("Microsatellite Instability-High");

    @NotNull
    public Map<Feature, String> extractSignatures(@NotNull ViccEntry viccEntry) {
        Map<Feature, String> signaturesPerFeature = Maps.newHashMap();
        if (viccEntry.source() == ViccSource.ONCOKB) {
            for (Feature feature : viccEntry.features()) {
                if (ONCOKB_SIGNATURES.contains(feature.name())) {
                    signaturesPerFeature.put(feature, feature.name());
                }
            }
        }
        return signaturesPerFeature;
    }

    @NotNull
    public static Signatures determineSignatures(@NotNull ViccSource source, @NotNull String typeEvent) {
        return ImmutableSignatures.builder().eventType(typeEvent).source(source.displayString()).sourceLink(source.toString()).build();
    }
}
