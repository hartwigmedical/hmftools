package com.hartwig.hmftools.ckb.classification;

import com.hartwig.hmftools.common.serve.classification.EventPreprocessor;

import org.jetbrains.annotations.NotNull;

public class ProteinAnnotationExtractor implements EventPreprocessor {

    public ProteinAnnotationExtractor() {
    }

    @NotNull
    @Override
    public String apply(@NotNull String event) {
        return event;
    }
}
