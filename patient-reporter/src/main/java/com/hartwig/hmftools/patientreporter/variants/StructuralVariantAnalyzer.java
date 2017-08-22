package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.svannotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.svannotation.StructuralVariantAnnotator;

import org.jetbrains.annotations.NotNull;

public class StructuralVariantAnalyzer {

    private final StructuralVariantAnnotator annotator;

    public StructuralVariantAnalyzer(final StructuralVariantAnnotator annotator) {
        this.annotator = annotator;
    }

    @NotNull
    public StructuralVariantAnalysis run(@NotNull final List<StructuralVariant> variants) {
        final List<StructuralVariantAnnotation> annotations = annotator.annotateVariants(variants);
        return new StructuralVariantAnalysis(annotations);
    }
}
