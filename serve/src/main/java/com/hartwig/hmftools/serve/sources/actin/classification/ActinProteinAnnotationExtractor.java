package com.hartwig.hmftools.serve.sources.actin.classification;

import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.serve.classification.EventPreprocessor;

import org.jetbrains.annotations.NotNull;

public class ActinProteinAnnotationExtractor implements EventPreprocessor {

    public ActinProteinAnnotationExtractor() {
    }

    @NotNull
    @Override
    public String apply(@NotNull String event) {
        return AminoAcids.forceSingleLetterProteinAnnotation(event);
    }
}

