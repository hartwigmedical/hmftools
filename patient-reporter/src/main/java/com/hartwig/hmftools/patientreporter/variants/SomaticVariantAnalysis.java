package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class SomaticVariantAnalysis {

    @NotNull
    public abstract List<EnrichedSomaticVariant> variantsToReport();

    public abstract double indelsPerMb();

    public abstract int mutationalLoad();

    public abstract double tumorMutationalBurden();
}
