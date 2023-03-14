package com.hartwig.hmftools.datamodel.purple;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = {NotNull.class, Nullable.class})
public abstract class PurpleCopyNumber {

    @NotNull
    public abstract String chromosome();

    public abstract int start();

    public abstract int end();

    public abstract double averageTumorCopyNumber();
}