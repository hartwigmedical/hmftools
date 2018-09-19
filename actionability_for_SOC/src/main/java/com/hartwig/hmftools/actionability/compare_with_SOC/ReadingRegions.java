package com.hartwig.hmftools.actionability.compare_with_SOC;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })

public abstract class ReadingRegions {

    @NotNull
    abstract String gene();

    @NotNull
    abstract String refSeq();

    @NotNull
    abstract String region();

    @NotNull
    abstract String amplificationDetection();
}
