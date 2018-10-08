package com.hartwig.hmftools.actionability.soc;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class ReadingPlatform {

    @NotNull
    abstract String chromosome();

    @NotNull
    abstract Integer startPosition();

    @NotNull
    abstract Integer endPosition();

    @NotNull
    abstract String targetRegionID();

    @NotNull
    abstract String target();
}
