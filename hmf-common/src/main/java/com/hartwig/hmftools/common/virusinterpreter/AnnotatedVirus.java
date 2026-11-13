package com.hartwig.hmftools.common.virusinterpreter;

import com.hartwig.hmftools.common.viridian.VirusType;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakendQCStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class AnnotatedVirus
{
    public abstract int taxid();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract VirusBreakendQCStatus qcStatus();

    public abstract int integrations();

    @Nullable
    public abstract VirusType interpretation();

    public abstract double percentageCovered();

    public abstract double meanCoverage();

    @Nullable
    public abstract Double expectedClonalCoverage();

    public abstract boolean reported();

    // Nullable for backwards compatibility
    @Nullable
    public abstract Boolean blacklisted();

    @NotNull
    public abstract VirusLikelihoodType virusDriverLikelihoodType();
}