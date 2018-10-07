package com.hartwig.hmftools.actionability.compare_with_SOC;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class ReadingRegionsBedFile {

    @NotNull
    abstract String chromosome();

    @NotNull
    abstract String startPosition();

    @NotNull
    abstract String endposition();

    @NotNull
    abstract String gene();

    @NotNull
    abstract String strand();
}
