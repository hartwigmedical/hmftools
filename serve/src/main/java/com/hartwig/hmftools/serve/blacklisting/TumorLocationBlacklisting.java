package com.hartwig.hmftools.serve.blacklisting;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class TumorLocationBlacklisting {

    @Nullable
    public abstract String blacklistCancerType();

    @Nullable
    public abstract String blacklistedDoid();
}
