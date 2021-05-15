package com.hartwig.hmftools.patientreporter.virusbreakend;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class VirusSummaryEntry {

    public abstract int taxidSpeciesId();

    @NotNull
    public abstract String reportVirusSummaryName();
}
