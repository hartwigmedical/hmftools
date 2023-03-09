package com.hartwig.hmftools.datamodel.purple;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = {NotNull.class, Nullable.class})
public abstract class PurpleGeneCopyNumber {

    @NotNull
    public abstract String chromosome();

    @NotNull
    public abstract String chromosomeBand();

    @NotNull
    public abstract String gene();

    @Nullable
    public abstract Double minCopyNumber();

    @Nullable
    public abstract Double minMinorAlleleCopyNumber();

}
