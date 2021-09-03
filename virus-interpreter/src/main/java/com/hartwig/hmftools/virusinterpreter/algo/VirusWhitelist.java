package com.hartwig.hmftools.virusinterpreter.algo;

import com.hartwig.hmftools.common.virus.VirusInterpretation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class VirusWhitelist {

    public abstract Integer taxidSpecies();

    public abstract boolean reportOnSummary();

    @NotNull
    public abstract VirusInterpretation virusInterpretation();

    @NotNull
    public abstract String nameSpecies();

    @Nullable
    public abstract Integer integratedMinimalCoverage();

    @Nullable
    public abstract Integer nonintegratedMinimalCoverage();

    @Nullable
    public abstract String integratedMeanDepth();

    @Nullable
    public abstract String nonintegratedMeanDepth();


}
