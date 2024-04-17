package com.hartwig.hmftools.datamodel.purple;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleGeneCopyNumber
{
    @NotNull
    String gene();

    @NotNull
    String chromosome();

    @NotNull
    String chromosomeBand();

    double minCopyNumber();

    double maxCopyNumber();

    double minMinorAlleleCopyNumber();
}
