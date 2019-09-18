package com.hartwig.hmftools.patientreporter;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class SampleMetadata {

    @NotNull
    public abstract String refSampleId();

    @NotNull
    public abstract String refSampleBarcode();

    @NotNull
    public abstract String tumorSampleId();

    @NotNull
    public abstract String tumorSampleBarcode();
}
