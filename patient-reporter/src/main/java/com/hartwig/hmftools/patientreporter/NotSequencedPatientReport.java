package com.hartwig.hmftools.patientreporter;

import java.util.Optional;

import com.hartwig.hmftools.patientreporter.algo.NotSequenceableReason;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableStudy;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class NotSequencedPatientReport implements PatientReport {
    @Override
    @NotNull
    public abstract SampleReport sampleReport();

    @NotNull
    public abstract NotSequenceableReason reason();

    @NotNull
    public abstract NotSequenceableStudy study();

    @Override
    @NotNull
    public abstract Optional<String> comments();
}

