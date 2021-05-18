package com.hartwig.hmftools.patientreporter.virusbreakend;


import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableVirusBreakend {

    @NotNull
    public abstract String virusName();

    public abstract int integrations();

    @Nullable
    public abstract String interpretations();
}
