package com.hartwig.hmftools.svannotation;

import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

public class VariantAnnotation {

    private final StructuralVariant Variant;
    private Breakend StartAnnotations;
    private Breakend EndAnnotations;

    public VariantAnnotation(final StructuralVariant variant) {
        Variant = variant;
        StartAnnotations = new Breakend(this);
        EndAnnotations = new Breakend(this);
    }

    public void setBreakendAnnotations(final Breakend start, final Breakend end) {
        StartAnnotations = start;
        EndAnnotations = end;
    }

    public StructuralVariant getVariant() {
        return Variant;
    }

    public Breakend getStart() {
        return StartAnnotations;
    }

    public Breakend getEnd() {
        return EndAnnotations;
    }

}
