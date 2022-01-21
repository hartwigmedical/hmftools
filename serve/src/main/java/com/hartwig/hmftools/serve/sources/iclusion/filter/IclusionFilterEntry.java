package com.hartwig.hmftools.serve.sources.iclusion.filter;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class IclusionFilterEntry {

    @NotNull
    public abstract IclusionFilterType type();

    @NotNull
    public abstract String value();
}
