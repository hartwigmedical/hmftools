package com.hartwig.hmftools.common.actionability.cnv;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionabilityCNVs {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String cnvType();

    @NotNull
    public abstract String source();

    @NotNull
    public abstract String reference();

    @NotNull
    public abstract String drugsName();

    @NotNull
    public abstract String drugsType();

    @NotNull
    public abstract String cancerType();

    @NotNull
    public abstract String hmfLevel();

    @NotNull
    public abstract String hmfResponse();
}
