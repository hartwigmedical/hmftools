package com.hartwig.hmftools.vicc.datamodel.jaxtrials;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class JaxTrialsIndication {

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String source();

    @NotNull
    public abstract String id();
}
