package com.hartwig.hmftools.common.linx;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableGeneDisruption {

    @NotNull
    public abstract String location();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String transcriptId();

    public abstract boolean canonical();

    @NotNull
    public abstract String range();

    @NotNull
    public abstract String type();

    @Nullable
    public abstract Double junctionCopyNumber();

    public abstract double undisruptedCopyNumber();

    public abstract int firstAffectedExon();

    public abstract int svId();

    @Nullable
    public abstract Integer clusterId();
}
