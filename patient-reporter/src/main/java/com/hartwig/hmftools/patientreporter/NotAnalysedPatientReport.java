package com.hartwig.hmftools.patientreporter;

import java.util.Optional;

import com.hartwig.hmftools.patientreporter.qcfail.NotAnalysableReason;
import com.hartwig.hmftools.patientreporter.qcfail.NotAnalysableStudy;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class NotAnalysedPatientReport implements PatientReport {
    @Override
    @NotNull
    public abstract SampleReport sampleReport();

    @NotNull
    public abstract NotAnalysableReason reason();

    @NotNull
    public abstract NotAnalysableStudy study();

    @Override
    @NotNull
    public abstract Optional<String> comments();

    @NotNull
    @Override
    public abstract String signaturePath();

    @NotNull
    @Override
    public abstract String logoRVA();
}

