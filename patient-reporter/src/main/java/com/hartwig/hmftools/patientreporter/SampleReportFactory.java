package com.hartwig.hmftools.patientreporter;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.hospital.HospitalModel;
import com.hartwig.hmftools.common.lims.Lims;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SampleReportFactory {

    private SampleReportFactory() {
    }

    @NotNull
    public static SampleReport fromLimsAndHospitalModel(@NotNull String tumorSample, @Nullable String refSample, @NotNull Lims lims,
            @NotNull HospitalModel hospitalModel, @Nullable PatientTumorLocation patientTumorLocation) {
        ImmutableSampleReport.Builder builder = ImmutableSampleReport.builder();

        // Ref sample is not resolved and also not relevant when generating QC Fail reports.
        if (refSample != null) {
            builder.refArrivalDate(lims.arrivalDate(refSample));
        }

        return builder.sampleId(tumorSample)
                .patientTumorLocation(patientTumorLocation)
                .refBarcode(lims.refBarcode(tumorSample))
                .tumorBarcode(lims.tumorBarcode(tumorSample))
                .tumorArrivalDate(lims.arrivalDate(tumorSample))
                .purityShallowSeq(lims.purityShallowSeq(tumorSample))
                .pathologyTumorPercentage(lims.pathologyTumorPercentage(tumorSample))
                .labProcedures(lims.labProcedures(tumorSample))
                .addressee(hospitalModel.addresseeStringForSample(tumorSample, lims.requesterName(tumorSample)))
                .projectName(lims.projectName(tumorSample))
                .requesterName(lims.requesterName(tumorSample))
                .requesterEmail(lims.requesterEmail(tumorSample))
                .submissionId(lims.submissionId(tumorSample))
                .hospitalPatientId(lims.hospitalPatientId(tumorSample))
                .hospitalPathologySampleId(lims.hospitalPathologySampleId(tumorSample))
                .build();
    }
}
