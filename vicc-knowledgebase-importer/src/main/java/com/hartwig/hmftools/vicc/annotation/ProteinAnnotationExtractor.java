package com.hartwig.hmftools.vicc.annotation;

import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.jetbrains.annotations.NotNull;

public final class ProteinAnnotationExtractor {

    private ProteinAnnotationExtractor() {
    }

    @NotNull
    public static String proteinAnnotation(@NotNull Feature feature) {
        return HotspotClassifier.extractProteinAnnotation(feature.name());
    }
}
