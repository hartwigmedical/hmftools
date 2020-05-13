package com.hartwig.hmftools.common.hospital;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class HospitalData {

    @NotNull
    public abstract String hospitalId();

    @NotNull
    public abstract String hospitalPI();

    @Nullable
    public abstract String requestName();

    @Nullable
    public abstract String requestEmail();

}
