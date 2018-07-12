package com.hartwig.hmftools.breakpointinspector.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class Range {

    public abstract int start();

    public abstract int end();

    @NotNull
    public static Range invert(@NotNull Range range) {
        return ImmutableRange.of(-range.end(), range.start());
    }
}
