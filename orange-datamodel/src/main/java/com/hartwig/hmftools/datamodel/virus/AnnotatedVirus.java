package com.hartwig.hmftools.datamodel.virus;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
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
    public abstract VirusInterpretation interpretation();

    public abstract double percentageCovered();

    public abstract double meanCoverage();

    @Nullable
    public abstract Double expectedClonalCoverage();

    public abstract boolean reported();

    @NotNull
    public abstract VirusLikelihoodType virusDriverLikelihoodType();
}