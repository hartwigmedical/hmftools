package com.hartwig.hmftools.svannotation;

import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

public class StructuralVariantAnnotation {

    private final StructuralVariant Variant;
    private BreakendAnnotations StartAnnotations;
    private BreakendAnnotations EndAnnotations;

    StructuralVariantAnnotation(final StructuralVariant variant) {
        Variant = variant;
        StartAnnotations = new BreakendAnnotations(this);
        EndAnnotations = new BreakendAnnotations(this);
    }

    void setBreakendAnnotations(final BreakendAnnotations start, final BreakendAnnotations end) {
        StartAnnotations = start;
        EndAnnotations = end;
    }

    public StructuralVariant getVariant() {
        return Variant;
    }

    public BreakendAnnotations getStartAnnotations() {
        return StartAnnotations;
    }

    public BreakendAnnotations getEndAnnotations() {
        return EndAnnotations;
    }
}
