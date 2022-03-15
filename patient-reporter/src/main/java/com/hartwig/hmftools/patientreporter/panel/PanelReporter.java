package com.hartwig.hmftools.patientreporter.panel;

import java.io.IOException;
import java.util.Optional;

import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumorFunctions;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.pipeline.PipelineVersionFile;
import com.hartwig.hmftools.patientreporter.QsFormNumber;
import com.hartwig.hmftools.patientreporter.SampleMetadata;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.SampleReportFactory;
import com.hartwig.hmftools.patientreporter.pipeline.PipelineVersion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PanelReporter {

    private static final Logger LOGGER = LogManager.getLogger(PanelReporter.class);

    @NotNull
    private final QCFailPanelReportData reportData;
    private final String reportDate;

    public PanelReporter(@NotNull final QCFailPanelReportData reportData, @NotNull final String reportDate) {
        this.reportData = reportData;
        this.reportDate = reportDate;
    }

    @NotNull
    public PanelReport run(@NotNull SampleMetadata sampleMetadata, @Nullable String comments, boolean correctedReport,
            boolean correctedReportExtern, @NotNull String expectedPipelineVersion, boolean overridePipelineVersion,
            @Nullable String pipelineVersionFile, boolean requirePipelineVersionFile, @NotNull String panelVCFname,
            @NotNull String panelGbase, @NotNull String panelQ30) throws IOException {

        String patientId = reportData.limsModel().patientId(sampleMetadata.tumorSampleBarcode());

        PatientPrimaryTumor patientPrimaryTumor =
                PatientPrimaryTumorFunctions.findPrimaryTumorForPatient(reportData.patientPrimaryTumors(), patientId);
        SampleReport sampleReport = SampleReportFactory.fromLimsModel(sampleMetadata, reportData.limsModel(), patientPrimaryTumor);

        String pipelineVersion = null;
        if (requirePipelineVersionFile) {
            assert pipelineVersionFile != null;
            pipelineVersion = PipelineVersionFile.majorDotMinorVersion(pipelineVersionFile);
            PipelineVersion.checkPipelineVersion(pipelineVersion, expectedPipelineVersion, overridePipelineVersion);
        }

        LimsCohortConfig cohort = sampleReport.cohort();

        if (cohort.cohortId().isEmpty()) {
            throw new IllegalStateException("QC fail report not supported for non-cancer study samples: " + sampleMetadata.tumorSampleId());
        }

        return ImmutablePanelReport.builder()
                .sampleReport(sampleReport)
                .qsFormNumber(QsFormNumber.FOR_080.display())
                .pipelineVersion(pipelineVersion)
                .VCFFilename(panelVCFname)
                .sampleGbase(panelGbase)
                .sampleQ30Value(panelQ30)
                .comments(Optional.ofNullable(comments))
                .isCorrectedReport(correctedReport)
                .isCorrectedReportExtern(correctedReportExtern)
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .reportDate(reportDate)
                .build();
    }
}