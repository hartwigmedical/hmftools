package com.hartwig.hmftools.patientreporter.qcfail;

import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.patientreporter.BaseReportData;
import com.hartwig.hmftools.patientreporter.ImmutableQCFailReport;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.hartwig.hmftools.patientreporter.SampleReport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class QCFailReporter {

    @NotNull
    abstract BaseReportData baseReportData();

    public QCFailReport run(@NotNull final String sample, @NotNull final QCFailReason reason, @Nullable final String comments) {
        final Lims lims = baseReportData().limsModel();

        final QCFailStudy study = QCFailStudy.fromSample(sample);

        assert study != null;

        final PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(baseReportData().patientTumorLocations(), sample);

        final SampleReport sampleReport = ImmutableSampleReport.builder()
                .sampleId(sample)
                .barcodeTumor(lims.tumorBarcode(sample))
                .barcodeReference(lims.refBarcode(sample))
                .patientTumorLocation(patientTumorLocation)
                .purityShallowSeq(lims.purityShallowSeq(sample))
                .pathologyTumorPercentage(lims.pathologyTumorPercentage(sample))
                .tumorArrivalDate(lims.arrivalDate(sample))
                .labProcedures(lims.labProcedures(sample))
                .addressee(baseReportData().hospitalModel().addresseeStringForSample(sample, lims.requesterName(sample)))
                .projectName(lims.projectName(sample))
                .requesterName(lims.requesterName(sample))
                .requesterEmail(lims.requesterEmail(sample))
                .submissionId(lims.submissionId(sample))
                .hospitalPatientId(lims.hospitalPatientId(sample))
                .hospitalPathologySampleId(lims.hospitalPathologySampleId(sample))
                .build();

        return ImmutableQCFailReport.of(sampleReport,
                reason,
                study,
                Optional.ofNullable(comments),
                baseReportData().signaturePath(),
                baseReportData().logoRVAPath());
    }
}
