package com.hartwig.hmftools.svannotation.analysis;

import java.util.List;

import com.hartwig.hmftools.svannotation.annotations.GeneDisruption;
import com.hartwig.hmftools.svannotation.annotations.GeneFusion;
import com.hartwig.hmftools.svannotation.annotations.StructuralVariantAnnotation;

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
