package com.hartwig.hmftools.patientreporter.qcfail;

import java.util.List;
import java.util.Optional;
import java.util.Set;

import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
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

    @Override
    @NotNull
    public abstract String qsFormNumber();

    @NotNull
    public abstract QCFailReason reason();

    @Nullable
    public abstract String wgsPurityString();

    @Nullable
    public abstract Set<PurpleQCStatus> purpleQC();

    @Override
    @NotNull
    public abstract Optional<String> comments();

    @Override
    public abstract boolean isCorrectedReport();

    @Override
    public abstract boolean isCorrectedReportExtern();

    @Override
    @NotNull
    public abstract String signaturePath();

    @Override
    @NotNull
    public abstract String logoRVAPath();

    @Override
    @NotNull
    public abstract String logoCompanyPath();

    @NotNull
    @Override
    public abstract List<PeachGenotype> peachGenotypes();

    @NotNull
    @Override
    public abstract String reportDate();
}