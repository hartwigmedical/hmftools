package com.hartwig.hmftools.iclusion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class IclusionCredentials {

    @NotNull
    public abstract String clientId();

    @NotNull
    public abstract String clientSecret();

    @NotNull
    public abstract String username();

    @NotNull
    public abstract String password();
}
