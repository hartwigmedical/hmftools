package com.hartwig.hmftools.common.virus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class InterpretedVirus {

    public abstract int taxid();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract VirusBreakendQCStatus qcStatus();

    public abstract int integrations();

    @Nullable
    public abstract VirusInterpretation interpretation();

    public abstract boolean reported();
}
