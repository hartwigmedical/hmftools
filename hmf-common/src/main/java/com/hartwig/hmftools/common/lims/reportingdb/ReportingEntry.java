package com.hartwig.hmftools.common.lims.reportingdb;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportingEntry {

    @NotNull
    public abstract String tumorBarcode();

    @NotNull
    public abstract String sampleId();

    @NotNull
    public abstract String cohort();

    @NotNull
    public abstract String reportDate();

    @NotNull
    public abstract String reportType();

    @NotNull
    public abstract String purity();

    @NotNull
    public abstract String hasReliableQuality();

    @NotNull
    public abstract String hasReliablePurity();
}
