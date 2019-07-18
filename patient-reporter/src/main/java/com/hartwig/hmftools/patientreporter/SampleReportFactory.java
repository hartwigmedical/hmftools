package com.hartwig.hmftools.patientreporter;

import java.time.LocalDate;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.hospital.HospitalModel;
import com.hartwig.hmftools.common.hospital.HospitalQuery;
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
    public static SampleReport fromLimsAndHospitalModel(@NotNull String tumorSample, @NotNull String refSample, @NotNull Lims lims,
            @NotNull HospitalModel hospitalModel, @Nullable PatientTumorLocation patientTumorLocation) {
        LocalDate arrivalDateRefSample = lims.arrivalDate(refSample);
        if (arrivalDateRefSample == null) {
            LOGGER.warn("Could not find arrival date for ref sample: " + refSample);
        }

        LocalDate arrivalDateTumorSample = lims.arrivalDate(tumorSample);
        if (arrivalDateTumorSample == null) {
            LOGGER.warn("Could not find arrival date for tumor sample: " + tumorSample);
        }

        HospitalQuery hospitalQuery = hospitalModel.queryHospitalDataForSample(tumorSample);

        return ImmutableSampleReport.builder().sampleId(tumorSample)
                .patientTumorLocation(patientTumorLocation)
                .refBarcode(lims.refBarcode(tumorSample))
                .refArrivalDate(arrivalDateRefSample)
                .tumorBarcode(lims.tumorBarcode(tumorSample))
                .tumorArrivalDate(arrivalDateTumorSample)
                .purityShallowSeq(lims.purityShallowSeq(tumorSample))
                .pathologyTumorPercentage(lims.pathologyTumorPercentage(tumorSample))
                .labProcedures(lims.labProcedures(tumorSample))
                .addressee(hospitalQuery.fullAddresseeString())
                .hospitalName(hospitalQuery.hospitalName())
                .hospitalPIName(hospitalQuery.principalInvestigatorName())
                .hospitalPIEmail(hospitalQuery.principalInvestigatorEmail())
                .projectName(lims.projectName(tumorSample))
                .requesterName(lims.requesterName(tumorSample))
                .requesterEmail(lims.requesterEmail(tumorSample))
                .submissionId(lims.submissionId(tumorSample))
                .hospitalPatientId(lims.hospitalPatientId(tumorSample))
                .hospitalPathologySampleId(lims.hospitalPathologySampleId(tumorSample))
                .build();
    }
}
