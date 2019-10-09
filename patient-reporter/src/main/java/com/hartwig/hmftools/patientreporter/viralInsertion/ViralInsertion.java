package com.hartwig.hmftools.patientreporter.viralInsertion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ViralInsertion {

    @NotNull
    public abstract String virus();

    @NotNull
    public abstract String countVirus();
}
