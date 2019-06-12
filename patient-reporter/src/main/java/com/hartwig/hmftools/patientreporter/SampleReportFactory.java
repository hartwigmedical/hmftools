package com.hartwig.hmftools.patientreporter;

import java.time.LocalDate;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.hospital.HospitalModel;
import com.hartwig.hmftools.common.lims.Lims;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SampleReportFactory {

    private static final Logger LOGGER = LogManager.getLogger(SampleReportFactory.class);

    private SampleReportFactory() {
    }

    @NotNull
    public static SampleReport fromLimsAndHospitalModel(@NotNull String tumorSample, @Nullable String refSample, @NotNull Lims lims,
            @NotNull HospitalModel hospitalModel, @Nullable PatientTumorLocation patientTumorLocation) {
        ImmutableSampleReport.Builder builder = ImmutableSampleReport.builder();

        // Ref sample is not resolved and also not relevant when generating QC Fail reports.
        if (refSample != null) {
            LocalDate arrivalDateRefSample = lims.arrivalDate(refSample);
            if (arrivalDateRefSample == null) {
                LOGGER.warn("Could not find arrival date for ref sample: " + refSample);
            }
            builder.refArrivalDate(arrivalDateRefSample);
        }

        LocalDate arrivalDateTumorSample = lims.arrivalDate(tumorSample);
        if (arrivalDateTumorSample == null) {
            LOGGER.warn("Could not find arrival date for tumor sample: " + tumorSample);
        }

        return builder.sampleId(tumorSample)
                .patientTumorLocation(patientTumorLocation)
                .refBarcode(lims.refBarcode(tumorSample))
                .tumorBarcode(lims.tumorBarcode(tumorSample))
                .tumorArrivalDate(arrivalDateTumorSample)
                .purityShallowSeq(lims.purityShallowSeq(tumorSample))
                .pathologyTumorPercentage(lims.pathologyTumorPercentage(tumorSample))
                .labProcedures(lims.labProcedures(tumorSample))
                .addressee(hospitalModel.fullAddresseeString(tumorSample))
                .hospitalName(hospitalModel.externalHospitalName(tumorSample))
                .hospitalPIName(hospitalModel.PIName(tumorSample))
                .hospitalPIEmail(hospitalModel.PIEmail(tumorSample))
                .projectName(lims.projectName(tumorSample))
                .requesterName(lims.requesterName(tumorSample))
                .requesterEmail(lims.requesterEmail(tumorSample))
                .submissionId(lims.submissionId(tumorSample))
                .hospitalPatientId(lims.hospitalPatientId(tumorSample))
                .hospitalPathologySampleId(lims.hospitalPathologySampleId(tumorSample))
                .build();
    }
}
