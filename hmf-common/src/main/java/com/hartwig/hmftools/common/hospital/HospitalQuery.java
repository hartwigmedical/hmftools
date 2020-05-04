package com.hartwig.hmftools.common.hospital;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class HospitalQuery {

    @NotNull
    public abstract String hospitalPA();

    @NotNull
    public abstract String analyseRequestName();

    @NotNull
    public abstract String analyseRequestEmail();

    @NotNull
    public abstract String hospitalId();

    @NotNull
    public abstract String hospitalName();

    @NotNull
    public abstract String hospitalAdres();
}
