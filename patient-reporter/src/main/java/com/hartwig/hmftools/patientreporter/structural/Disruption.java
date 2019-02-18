package com.hartwig.hmftools.patientreporter.structural;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class Disruption {

    public abstract boolean reportable();

    @NotNull
    public abstract String svId();

    @NotNull
    public abstract String chromosome();

    @NotNull
    public abstract String position();

    public abstract int orientation();

    @NotNull
    public abstract String type();

    @Nullable
    public abstract Double ploidy();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String chrBand();

    @NotNull
    public abstract String transcript();

    public abstract int strand();

    @NotNull
    public abstract String regionType();

    @NotNull
    public abstract String codingType();

    public abstract boolean canonical();

    @NotNull
    public abstract String biotype();

    public abstract int exonUp();

    public abstract int exonDown();

    public abstract boolean isDisruptive();
}
