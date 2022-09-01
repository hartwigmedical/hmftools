package com.hartwig.hmftools.patientreporter.algo.orange.cuppa;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CuppaPrediction {

    @NotNull
    public abstract String cancerType();

    public abstract double likelihood();
}

