package com.hartwig.hmftools.vicc.datamodel.jax;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class JaxMolecularProfile {

    @NotNull
    public abstract String profileName();

    @NotNull
    public abstract String id();
}
