package com.hartwig.hmftools.datamodel.purple;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = {NotNull.class, Nullable.class})
public abstract class PurpleGeneCopyNumber {

    @NotNull
    public abstract String chromosome();

    @NotNull
    public abstract String chromosomeBand();

    @NotNull
    public abstract String geneName();

    @Nullable
    public abstract Double minCopyNumber();

    @Nullable
    public abstract Double minMinorAlleleCopyNumber();

}
