package com.hartwig.hmftools.common.lims.hospital;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class HospitalContact {

    @NotNull
    public abstract String hospitalId();

    @NotNull
    public abstract String hospitalPI();

    @Nullable
    public abstract String requesterName();

    @Nullable
    public abstract String requesterEmail();

}
