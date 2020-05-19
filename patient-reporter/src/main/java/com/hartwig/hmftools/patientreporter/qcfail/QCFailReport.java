package com.hartwig.hmftools.patientreporter.qcfail;

import java.util.Optional;

import com.hartwig.hmftools.common.lims.LimsStudy;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class QCFailReport implements PatientReport {

    @Override
    @NotNull
    public abstract SampleReport sampleReport();

    @NotNull
    public abstract QCFailReason reason();

    @NotNull
    public abstract LimsStudy study();

    @Nullable
    public abstract String wgsPurityString();

    @Override
    @NotNull
    public abstract Optional<String> comments();

    @Override
    public abstract boolean isCorrectedReport();

    @Override
    @NotNull
    public abstract String signaturePath();

    @Override
    @NotNull
    public abstract String logoRVAPath();

    @Override
    @NotNull
    public abstract String logoCompanyPath();
}

