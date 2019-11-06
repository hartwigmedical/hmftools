package com.hartwig.hmftools.common.reportablegenomicalterations;

import com.hartwig.hmftools.common.bachelor.GermlineVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableGermlineVariant {

    @NotNull
    public abstract GermlineVariant variant();

    public abstract double driverLikelihood();
}
