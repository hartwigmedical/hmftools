package com.hartwig.hmftools.patientreporter;

import java.util.Optional;

import org.jetbrains.annotations.NotNull;

public interface PatientReport {

    @NotNull
    SampleReport sampleReport();

    @NotNull
    default String user() {
        return System.getProperty("user.name");
    }

    @NotNull
    Optional<String> comments();

    @NotNull
    String signaturePath();

    @NotNull
    String logoRVA();
}
