package com.hartwig.hmftools.patientreporter.reportingdb;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportingEntry {

    @NotNull
    public abstract String sampleId();

    @NotNull
    public abstract String tumorBarcode();

    @NotNull
    public abstract String reportDate();

    @NotNull
    public abstract String sourceReport();

    @Nullable
    public abstract String purity();

    @Nullable
    public abstract String status();

    @Nullable
    public abstract String qcStatus();

}
