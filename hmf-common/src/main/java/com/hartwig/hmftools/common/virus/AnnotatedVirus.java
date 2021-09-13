package com.hartwig.hmftools.common.virus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class AnnotatedVirus {

    public abstract int taxid();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract VirusBreakendQCStatus qcStatus();

    public abstract int integrations();

    @NotNull
    public abstract String interpretation();

    public abstract double percentageCovered();

    public abstract double meanCoverage();

    @Nullable
    public abstract Double expectedClonalCoverage();

    public abstract boolean reported();

}
