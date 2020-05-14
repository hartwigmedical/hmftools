package com.hartwig.hmftools.common.hospital;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class HospitalData {

    @Nullable
    public abstract String hospitalPI();

    @Nullable
    public abstract String requesterName();

    @Nullable
    public abstract String requesterEmail();

    @Nullable
    public abstract String hospitalName();

    @Nullable
    public abstract String hospitalAddress();

}
