package com.hartwig.hmftools.patientreporter.algo;

import com.hartwig.hmftools.patientreporter.purple.PurpleAnalysis;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalysis;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class },
             allParameters = true)
abstract class GenomeAnalysis {
    @NotNull
    abstract String sample();

    @NotNull
    abstract VariantAnalysis variantAnalysis();

    @NotNull
    abstract PurpleAnalysis purpleAnalysis();

    @NotNull
    abstract StructuralVariantAnalysis structuralVariantAnalysis();

    abstract double indelsPerMb();
}
