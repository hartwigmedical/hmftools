package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.svannotation.StructuralVariantAnnotation;

import org.jetbrains.annotations.NotNull;

public class StructuralVariantAnalysis {

    @NotNull
    private List<StructuralVariantAnnotation> annotations;

    public StructuralVariantAnalysis(@NotNull final List<StructuralVariantAnnotation> annotations) {
        this.annotations = annotations;
    }

    @NotNull
    public List<StructuralVariantAnnotation> getAnnotations() {
        return ImmutableList.copyOf(annotations);
    }
}
