package com.hartwig.hmftools.common.variant.structural.annotation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
public abstract class ReportableDisruption
{
    @NotNull
    public abstract int svId();

    @NotNull
    public abstract String chromosome();

    @NotNull
    public abstract int orientation();

    @NotNull
    public abstract int strand();

    @NotNull
    public abstract String chrBand();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract boolean canonical();

    @NotNull
    public abstract String type();

    @Nullable
    public abstract Double ploidy();

    public abstract int exonUp();

    public abstract int exonDown();

}
