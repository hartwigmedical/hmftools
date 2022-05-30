package com.hartwig.hmftools.patientreporter;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class SampleMetadata {

    @Nullable
    public abstract String refSampleId();

    @Nullable
    public abstract String refSampleBarcode();

    @NotNull
    public abstract String tumorSampleId();

    @NotNull
    public abstract String tumorSampleBarcode();

    @NotNull
    public abstract String sampleNameForReport();
}
