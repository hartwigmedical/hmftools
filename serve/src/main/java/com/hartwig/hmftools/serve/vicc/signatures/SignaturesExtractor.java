package com.hartwig.hmftools.serve.vicc.signatures;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.vicc.ViccExtractor;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SignaturesExtractor {
    private static final Logger LOGGER = LogManager.getLogger(SignaturesExtractor.class);

    private static final Set<String> SIGNATURES =
            Sets.newHashSet("Microsatellite Instability-High");

    @NotNull
    public Map<Feature, String> extractSignatures(@NotNull ViccEntry viccEntry) {
        Map<Feature, String> signaturesPerFeature = Maps.newHashMap();
            for (Feature feature : viccEntry.features()) {
                if (SIGNATURES.contains(feature.name())) {
                    signaturesPerFeature.put(feature, feature.name());
                }
            }
        return signaturesPerFeature;
    }

    @NotNull
    public static Signatures determineSignatures(@NotNull ViccSource source, @NotNull String typeEvent) {
        return ImmutableSignatures.builder().eventType(typeEvent).source(source.displayString()).sourceLink(source.toString()).build();
    }
}
