package com.hartwig.hmftools.common.fusion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableDisruption
{
    public abstract int svId();

    @NotNull
    public abstract String chromosome();

    public abstract int orientation();

    public abstract int strand();

    @NotNull
    public abstract String chrBand();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String type();

    @Nullable
    public abstract Double ploidy();

    public abstract int exonUp();

    public abstract int exonDown();

    public abstract double undisruptedCopyNumber();

}
