package com.hartwig.hmftools.ckb.classification;

import com.hartwig.hmftools.common.serve.classification.EventPreprocessor;

import org.jetbrains.annotations.NotNull;

public class CkbProteinAnnotationExtractor implements EventPreprocessor {

    public CkbProteinAnnotationExtractor() {
    }

    @NotNull
    @Override
    public String apply(@NotNull String event) {
        String proteinAnnotation = event;
        // Cut out the trailing stop gained in case a stop gained is following on a frameshift
        int trailingStopGained = proteinAnnotation.indexOf("fs*");
        if (trailingStopGained > 0) {
            proteinAnnotation = proteinAnnotation.substring(0, trailingStopGained + 2);
        }

        return proteinAnnotation;
    }
}
