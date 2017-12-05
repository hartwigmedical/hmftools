package com.hartwig.hmftools.svannotation.annotations;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;

public class StructuralVariantAnnotation {

    @NotNull
    private final StructuralVariant Variant;
    private final List<GeneAnnotation> Annotations = Lists.newArrayList();

    public StructuralVariantAnnotation(@NotNull final StructuralVariant variant) {
        Variant = variant;
    }

    @NotNull
    public List<GeneAnnotation> getAnnotations() {
        return Annotations;
    }

    @NotNull
    public StructuralVariant getVariant() {
        return Variant;
    }

    @NotNull
    public List<GeneAnnotation> getStart() {
        return Annotations.stream().filter(GeneAnnotation::isStart).collect(Collectors.toList());
    }

    @NotNull
    public List<GeneAnnotation> getEnd() {
        return Annotations.stream().filter(a -> !a.isStart()).collect(Collectors.toList());
    }
}
