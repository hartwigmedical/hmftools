package com.hartwig.hmftools.patientreporter.qcfail;

import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.patientreporter.SampleMetadata;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.SampleReportFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class QCFailReporter {

    @NotNull
    private final QCFailReportData reportData;

    public QCFailReporter(@NotNull final QCFailReportData reportData) {
        this.reportData = reportData;
    }

    @NotNull
    public QCFailReport run(@NotNull SampleMetadata sampleMetadata, @NotNull QCFailReason reason, @Nullable String comments,
            boolean correctedReport) {
        QCFailStudy study = QCFailStudy.fromSampleId(sampleMetadata.tumorSampleId());

        assert study != null;

        PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(reportData.patientTumorLocations(),
                        sampleMetadata.tumorSampleId());

        SampleReport sampleReport = SampleReportFactory.fromLimsModel(sampleMetadata,
                reportData.limsModel(),
                patientTumorLocation);

        return ImmutableQCFailReport.builder()
                .sampleReport(sampleReport)
                .reason(reason)
                .study(study)
                .comments(Optional.ofNullable(comments))
                .isCorrectedReport(correctedReport)
                .signaturePath(reportData.signaturePath())
                .logoRVAPath(reportData.logoRVAPath())
                .logoCompanyPath(reportData.logoCompanyPath())
                .build();
    }
}
