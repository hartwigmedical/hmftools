package com.hartwig.hmftools.vicc.annotation;

import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ProteinAnnotationExtractor {

    private ProteinAnnotationExtractor() {
    }

    @NotNull
    public static String proteinAnnotation(@NotNull Feature feature) {
        String featureName = feature.name();
        return DetermineHotspot.isHotspot(featureName) ? DetermineHotspot.extractProteinAnnotation(featureName) : Strings.EMPTY;
    }
}
