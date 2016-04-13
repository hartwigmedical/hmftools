package com.hartwig.hmftools.boggs;

import org.jetbrains.annotations.NotNull;

public class PatientData {
    @NotNull
    private final SampleData refSample;
    @NotNull
    private final SampleData tumorSample;

    public PatientData(@NotNull SampleData refSample, @NotNull SampleData tumorSample) {
        this.refSample = refSample;
        this.tumorSample = tumorSample;
    }

    @Override
    public String toString() {
        return "PatientData{" +
                "refSample=" + refSample +
                ", tumorSample=" + tumorSample +
                '}';
    }
}
