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
    public static SampleReport fromLimsAndHospitalModel(@NotNull SampleMetadata sampleMetadata, @NotNull Lims lims,
            @NotNull HospitalModel hospitalModel, @Nullable PatientTumorLocation patientTumorLocation) {
        String refSampleId = sampleMetadata.refSampleId();
        String tumorSampleId = sampleMetadata.tumorSampleId();

        LocalDate arrivalDateRefSample = lims.arrivalDate(refSampleId);
        if (arrivalDateRefSample == null) {
            LOGGER.warn("Could not find arrival date for ref sample: " + refSampleId);
        }

        LocalDate arrivalDateTumorSample = lims.arrivalDate(tumorSampleId);
        if (arrivalDateTumorSample == null) {
            LOGGER.warn("Could not find arrival date for tumor sample: " + tumorSampleId);
        }

        HospitalQuery hospitalQuery = hospitalModel.queryHospitalDataForSample(tumorSampleId);

        return ImmutableSampleReport.builder()
                .sampleMetadata(sampleMetadata)
                .patientTumorLocation(patientTumorLocation)
                .refArrivalDate(arrivalDateRefSample)
                .tumorArrivalDate(arrivalDateTumorSample)
                .purityShallowSeq(lims.purityShallowSeq(tumorSampleId))
                .pathologyTumorPercentage(lims.pathologyTumorPercentage(tumorSampleId))
                .labProcedures(lims.labProcedures(tumorSampleId))
                .addressee(hospitalQuery.fullAddresseeString())
                .hospitalName(hospitalQuery.hospitalName())
                .hospitalPIName(hospitalQuery.principalInvestigatorName())
                .hospitalPIEmail(hospitalQuery.principalInvestigatorEmail())
                .projectName(lims.projectName(tumorSampleId))
                .requesterName(lims.requesterName(tumorSampleId))
                .requesterEmail(lims.requesterEmail(tumorSampleId))
                .submissionId(lims.submissionId(tumorSampleId))
                .hospitalPatientId(lims.hospitalPatientId(tumorSampleId))
                .hospitalPathologySampleId(lims.hospitalPathologySampleId(tumorSampleId))
                .build();
    }
}
