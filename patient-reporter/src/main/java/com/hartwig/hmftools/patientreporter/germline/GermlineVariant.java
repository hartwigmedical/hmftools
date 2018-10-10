package com.hartwig.hmftools.patientreporter.germline;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GermlineVariant {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String variant();

    @NotNull
    public abstract String impact();

    @NotNull
    public abstract String readDepth();

    @NotNull
    public abstract String ploidyVaf();

    @NotNull
    public abstract String germlineStatus();

    @NotNull
    public abstract String biallelic();

}
