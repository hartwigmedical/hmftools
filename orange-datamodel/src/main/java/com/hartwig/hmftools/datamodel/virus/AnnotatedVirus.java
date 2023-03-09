package com.hartwig.hmftools.datamodel.virus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class AnnotatedVirus {

    @NotNull
    public abstract String name();

    @NotNull
    public abstract VirusBreakendQCStatus qcStatus();

    public abstract int integrations();

    @Nullable
    public abstract String interpretation();

    public abstract double percentageCovered();

    public abstract boolean reported();
}