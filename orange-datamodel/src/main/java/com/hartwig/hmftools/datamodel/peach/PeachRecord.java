package com.hartwig.hmftools.datamodel.peach;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.Set;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = {NotNull.class, Nullable.class})
public abstract class PeachRecord {

    @NotNull
    public abstract Set<PeachGenotype> entries();
}
