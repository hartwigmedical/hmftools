package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class SampleData {

    @NotNull
    public abstract String sampleId();

    @NotNull
    public abstract LocalDate arrivalDate();

    @Nullable
    public abstract LocalDate samplingDate();

    @Nullable
    public abstract Integer dnaNanograms();

    @Nullable
    public abstract Double tumorPercentage();

    @NotNull
    @Value.Derived
    public LocalDate date() {
        final LocalDate samplingDate = samplingDate();
        return samplingDate != null ? samplingDate : arrivalDate();
    }

    @Override
    public String toString() {
        return String.format("Sample {%s}", sampleId());
    }
}
