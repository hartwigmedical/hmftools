package com.hartwig.hmftools.patientreporter.algo;

import java.util.Optional;

import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.common.cuppa.MolecularTissueOrigin;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class AnalysedPatientReport implements PatientReport {

    @Override
    @NotNull
    public abstract SampleReport sampleReport();

    @NotNull
    @Override
    public abstract String qsFormNumber();

    @NotNull
    public abstract String clinicalSummary();

    @Nullable
    public abstract String pipelineVersion();

    @NotNull
    public abstract GenomicAnalysis genomicAnalysis();

    @Nullable
    public abstract MolecularTissueOrigin molecularTissueOrigin();

    @NotNull
    public abstract String circosPath();

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
