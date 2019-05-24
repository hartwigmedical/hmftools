package com.hartwig.hmftools.patientreporter.qcfail;

import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.patientreporter.BaseReportData;
import com.hartwig.hmftools.patientreporter.ImmutableQCFailReport;
import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.SampleReportFactory;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class QCFailReporter {

    @NotNull
    abstract BaseReportData baseReportData();

    @NotNull
    public QCFailReport run(@NotNull final String sample, @NotNull final QCFailReason reason, @Nullable final String comments) {
        final QCFailStudy study = QCFailStudy.fromSample(sample);

        assert study != null;

        final PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(baseReportData().patientTumorLocations(), sample);

        final SampleReport sampleReport = SampleReportFactory.fromLimsAndHospitalModel(sample,
                null,
                baseReportData().limsModel(),
                baseReportData().hospitalModel(),
                patientTumorLocation);

        return ImmutableQCFailReport.of(sampleReport,
                reason,
                study,
                Optional.ofNullable(comments),
                baseReportData().signaturePath(),
                baseReportData().logoRVAPath(),
                baseReportData().logoCompanyPath());
    }
}
