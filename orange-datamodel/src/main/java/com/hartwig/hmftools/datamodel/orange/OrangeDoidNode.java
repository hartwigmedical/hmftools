package com.hartwig.hmftools.datamodel.orange;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeDoidNode
{
    @NotNull
    String doid();

    @Nullable
    String doidTerm();
}
