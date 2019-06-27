package com.hartwig.hmftools.patientreporter.qcfail;

import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.patientreporter.BaseReportData;
import com.hartwig.hmftools.patientreporter.ImmutableQCFailReport;
import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.SampleReportFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class QCFailReporter {

    @NotNull
    private final BaseReportData reportData;

    public QCFailReporter(@NotNull final BaseReportData reportData) {
        this.reportData = reportData;
    }

    @NotNull
    public QCFailReport run(@NotNull String sample, @NotNull QCFailReason reason, @Nullable String comments) {
        QCFailStudy study = QCFailStudy.fromSample(sample);

        assert study != null;

        PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(reportData.patientTumorLocations(), sample);

        SampleReport sampleReport = SampleReportFactory.fromLimsAndHospitalModel(sample,
                null,
                reportData.limsModel(),
                reportData.hospitalModel(),
                patientTumorLocation);

        return ImmutableQCFailReport.of(sampleReport,
                reason,
                study,
                Optional.ofNullable(comments),
                reportData.signaturePath(),
                reportData.logoRVAPath(),
                reportData.logoCompanyPath());
    }
}
