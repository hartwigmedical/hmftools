package com.hartwig.hmftools.common.variant.structural.annotation;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;

import org.jetbrains.annotations.NotNull;

public class StructuralVariantAnnotation {

    @NotNull
    private final StructuralVariantData variant;
    @NotNull
    private final List<GeneAnnotation> annotations = Lists.newArrayList();

    public StructuralVariantAnnotation(@NotNull final StructuralVariantData variant)
    {
        this.variant = variant;
    }

    @NotNull
    public StructuralVariantData variant() {
        return variant;
    }

    @NotNull
    public List<GeneAnnotation> annotations() {
        return annotations;
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
