package com.hartwig.hmftools.datamodel.isofox;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface NovelSpliceJunction
{
    @NotNull
    String gene();

    @NotNull
    String chromosome();

    int junctionStart();

    int junctionEnd();

    @NotNull
    AltSpliceJunctionType type();

    int fragmentCount();

    int depthStart();

    int depthEnd();

    @NotNull
    AltSpliceJunctionContext regionStart();

    @NotNull
    AltSpliceJunctionContext regionEnd();

    int cohortFrequency();
}
