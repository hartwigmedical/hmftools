package com.hartwig.hmftools.patientdb.readers.wide;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class WideBiopsyData {

    @NotNull
    public abstract String patientId();

    @NotNull
    public abstract String wideId();

    @NotNull
    public abstract String dataAvailable();

    @NotNull
    public abstract String tissueId();

    @NotNull
    public abstract String bioptDate();

    @NotNull
    public abstract String WGSsuccesfull();
}
