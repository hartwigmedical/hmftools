package com.hartwig.hmftools.common.hospital;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class HospitalCore {

    @NotNull
    public abstract String internalHospitalName();

    @NotNull
    public abstract String externalHospitalName();

    @NotNull
    public abstract String addressZip();

    @NotNull
    public abstract String addressCity();
}
