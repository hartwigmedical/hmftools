package com.hartwig.hmftools.serve.blacklisting;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class TumorLocationBlacklisting {

    @NotNull
    public abstract String blacklistCancerType();

    @NotNull
    public abstract String blacklistedDoid();
}
