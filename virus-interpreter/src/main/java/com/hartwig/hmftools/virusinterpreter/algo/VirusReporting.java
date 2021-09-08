package com.hartwig.hmftools.virusinterpreter.algo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class VirusReporting {

    public abstract boolean reportOnSummary();

    @NotNull
    public abstract String virusInterpretation();

    @Nullable
    public abstract Integer integratedMinimalCoverage();

    @Nullable
    public abstract Integer nonIntegratedMinimalCoverage();

}
