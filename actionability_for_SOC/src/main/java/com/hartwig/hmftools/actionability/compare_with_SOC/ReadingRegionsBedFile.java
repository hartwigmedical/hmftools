package com.hartwig.hmftools.actionability.compare_with_SOC;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })

public abstract class ReadingRegionsBedFile {

    @NotNull
    abstract String chromosome();

    @NotNull
    abstract Integer startPosition();

    @NotNull
    abstract Integer endposition();

    @NotNull
    abstract String gene();

    @NotNull
    abstract String strand();

}
