package com.hartwig.hmftools.common.variant.structural.annotation;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class StructuralVariantAnalysis {

    @NotNull
    public abstract List<StructuralVariantAnnotation> annotations();

    @NotNull
    public abstract List<GeneFusion> fusions();

    @NotNull
    public abstract List<GeneDisruption> disruptions();
}
