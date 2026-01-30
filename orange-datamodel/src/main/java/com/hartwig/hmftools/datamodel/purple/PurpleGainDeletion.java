package com.hartwig.hmftools.datamodel.purple;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleGainDeletion
{
    @NotNull
    PurpleDriver driver();

    @NotNull
    CopyNumberInterpretation interpretation();

    // only populated for germline amp dels
    @Nullable
    GermlineAmpDelFields germlineAmpDelFields();

    @NotNull
    default String gene()
    {
        return driver().gene();
    }

    @NotNull
    String chromosome();

    @NotNull
    String chromosomeBand();

    @NotNull
    default String transcript()
    {
        return driver().transcript();
    }

    default boolean isCanonical()
    {
        return driver().isCanonical();
    }

    // these are the tumor copies
    double minCopies();

    double maxCopies();

    double minMinorAlleleCopies();
}
