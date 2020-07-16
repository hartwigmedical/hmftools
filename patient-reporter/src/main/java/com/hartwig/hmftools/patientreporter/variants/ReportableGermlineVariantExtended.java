package com.hartwig.hmftools.patientreporter.variants;

import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableGermlineVariantExtended
{
    @NotNull
    public abstract ReportableGermlineVariant variant();

    public abstract double driverLikelihood();
}
