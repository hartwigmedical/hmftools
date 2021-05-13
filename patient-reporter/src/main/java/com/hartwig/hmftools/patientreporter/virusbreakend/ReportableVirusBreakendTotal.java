package com.hartwig.hmftools.patientreporter.virusbreakend;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableVirusBreakendTotal {

    @NotNull
    public abstract List<ReportableVirusBreakend> reportableViruses();

    @NotNull
    public abstract String virusNameSummary();
}
