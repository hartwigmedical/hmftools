package com.hartwig.hmftools.svannotation.annotations;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;

public class StructuralVariantAnnotation {

    @NotNull
    private final StructuralVariant variant;
    @NotNull
    private final List<GeneAnnotation> annotations = Lists.newArrayList();

    public StructuralVariantAnnotation(@NotNull final StructuralVariant variant) {
        this.variant = variant;
    }

    @NotNull
    public List<GeneAnnotation> annotations() {
        return annotations;
    }

    @NotNull
    public StructuralVariant variant() {
        return variant;
    }

    @NotNull
    public List<GeneAnnotation> start() {
        return annotations.stream().filter(GeneAnnotation::isStart).collect(Collectors.toList());
    }

    @NotNull
    public List<GeneAnnotation> end() {
        return annotations.stream().filter(annotation -> !annotation.isStart()).collect(Collectors.toList());
    }
}
