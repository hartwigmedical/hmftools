package com.hartwig.hmftools.patientreporter.purple;

import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableGainLoss {

    @NotNull
    public abstract CopyNumberInterpretation interpretation();

    @NotNull
    public abstract String chromosome();

    @NotNull
    public abstract String chromosomeBand();

    @NotNull
    public abstract String gene();

    public abstract long copies();
}
